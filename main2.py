import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
import scipy as sp
from scipy.io.wavfile import write

# Equation implemented from Chaigne and Askenfeld 1993
# Define constants

n_physical_segments = 50  # no. of string segments


# Values follow for C4 as listed in Chainge and Askenfelt 1993
fund_freq = 262  # Hz, assume standard tuning frequency A4
L = 0.62  # m
M_s = 3.93 * 3  # string mass
T = 670  # string tension
E = 2.0e11  # Youngs modulus
eps = 3.82e-5  # string stiffness parameter
b_1 = 0.5  # damping coefficient
b_3 = 6.25e-9  # damping coefficient

M_h = 2.97 * 3  # hammer mass
mu = M_h / L  # Transversal string desntiy, Placeholder value
# mu = 10e6
c = np.sqrt(T / mu)
v_h0 = 1.5  # m / s**2; note 4.0 1.5 0.5 for resp. forte mezzo forte piano
p = 2.5  # ideally between 2 and 3
K_h = 4.5e9  # generalised hammer stiffness
hammer_striking_position = 0.12
hammer_location_index = int(
    hammer_striking_position * n_physical_segments + 1
)  # add one due to the leftmost virtual point
index_to_phyiscal_position_factor = L / (n_physical_segments + 1) 
eta_zeroth_step = 0  # position of the hammer at t=0 random value?????????? 

dx = L / n_physical_segments
dt = 0.0001
max_t = 10.0


samp_freq = 32000  # Hz

D = 1 + b_1 * dt + 2 * b_3 / dt
r = c * dt / dx

a_1 = (2 - 2 * r**2 + b_3 / dt - 6 * eps * n_physical_segments**2 * r**2) / D
a_2 = (-1 + b_1 * dt + 2 * b_3 / dt) / D
a_3 = (r**2 * (1 + 4 * eps * n_physical_segments**2)) / D
a_4 = (b_3 / dt - eps * n_physical_segments**2 * r**2) / D
a_5 = (-b_3 / dt) / D


def hammer_force(eta, y_x0, K_h=K_h, p=p):
    """Note eta should be put in as a scalar"""
    # print("k_0 string extension as follows: " + str(y_x0))
    return K_h * np.abs(eta - y_x0) ** p


def update_virtual_points(array):
    """updates the virtual (ghost) points of an array.
    These points should be equal to the point directly
    neighbouring the attachment point

    Args:
        array (arr): input arrau
    """
    array[0] = -1 * array[2]
    array[-1] = -1 * array[-3]

def simulate_string():

    # we will need two virtual points to be able to perform the recursion
    n_virtual_points = 2
    n_steps = int(max_t / dt)
    n_total_segments = n_physical_segments + n_virtual_points

    string_position_arr = np.zeros([n_steps, n_total_segments])
    eta = np.zeros(n_steps)
    F_h = np.zeros(n_steps)


    # make a window function around the hammer location
    window_function = np.zeros(n_total_segments)
    window_function[hammer_location_index - 2 : hammer_location_index + 1] = [1.5, 0, 1.5]


    # initialize by doing 3 steps
    # first, make an estimate for the first step
    zeroth_step = np.zeros(n_total_segments)
    first_step = 0.5 * (np.roll(zeroth_step, 1) + np.roll(zeroth_step, -1))
    eta_first_step = v_h0 * dt
    force_magnitude_first_step = hammer_force(
        eta_first_step, first_step[hammer_location_index]
    )
    # calculate second step ignoring stiffness
    second_step = (
        np.roll(first_step, -1)
        + np.roll(first_step, 1)
        - zeroth_step
        + (dt * n_physical_segments * force_magnitude_first_step * window_function) / M_s
    )
    update_virtual_points(second_step)


    eta_second_step = (
        2 * eta_first_step - eta_zeroth_step - (dt**2 * force_magnitude_first_step)
    )
    force_magnitude_second_step = hammer_force(
        eta_second_step, second_step[hammer_location_index]
    )

    string_position_arr[0, :] = zeroth_step
    string_position_arr[1, :] = first_step
    string_position_arr[2, :] = second_step

    ############

    start_index = 3  # start after the initialization
    for i in range(start_index, n_steps - 1):
        two_steps_earlier = string_position_arr[i - 2, :]
        one_step_earlier = string_position_arr[i - 1, :]
        current = string_position_arr[i, :]

        next_step = (
            a_1 * current
            + a_2 * one_step_earlier
            + a_3 * (np.roll(current, 1) + np.roll(current, -1))
            + a_4 * (np.roll(current, 2) + np.roll(current, -2))
            + a_5
            * (
                np.roll(one_step_earlier, 1)
                + np.roll(one_step_earlier, -1)
                + two_steps_earlier
            )
            + 0  # force term
        )

        # the duplicate points mirror the points right before the boundry
        update_virtual_points(next_step)

        string_position_arr[i + 1, :] = next_step

    return string_position_arr

    
def make_fft_arr(string_position_arr, dt, max_t, string_position_index = 50):
    local_string_position = string_position_arr[:, string_position_index]
    samples_per_chunk = 100
    n_samples = len(local_string_position)
    n_samples_rounded_up = int(math.ceil(n_samples)) # round up for ease of operation
    splits = n_samples_rounded_up/samples_per_chunk

    n_steps = int(max_t / dt)
    freq = np.fft.fftfreq(samples_per_chunk, dt)
    output_arr = np.zeros((len(local_string_position), len(freq)))

    # split into arrays of 10 long except possibly the last split if not divisible
    chunks = np.array_split(local_string_position, splits)

    for i, chunk in enumerate(chunks):
        chunk_length = len(chunk)
        if chunk_length == samples_per_chunk:
            print(f'{samples_per_chunk=}')
            print(f'{chunk_length=}')
            print('right length')
            # for the regular chunks
            abs_fft = np.abs(np.fft.fft(chunk, axis=0))
            output_arr[i * samples_per_chunk: (i+1) * samples_per_chunk] = np.tile(abs_fft, (samples_per_chunk, 1))
        else:
            print(f'{samples_per_chunk=}')
            print(f'{chunk_length=}')
            print('wrong length')
            # for the regular chunks
            #for the last chunk thats possibly indivisble by samples_per_chunk
            output_arr[i * samples_per_chunk: i * samples_per_chunk + chunk_length] = 0

            
    positive_frequency_index = int(len(freq)/2) 
    
    return freq[:positive_frequency_index], output_arr[:,:positive_frequency_index]



    
def make_fft_arr2(string_position_arr, dt, max_t, string_position_index = 50):
    local_string_position = string_position_arr[:, string_position_index]
    samples_per_chunk = 1000
    n_samples = len(local_string_position)
    n_samples_rounded_up = int(math.ceil(n_samples)) # round up for ease of operation
    splits = n_samples_rounded_up/samples_per_chunk

    n_steps = int(max_t / dt)
    freq = np.fft.fftfreq(samples_per_chunk, dt)
    output_arr = np.zeros((len(local_string_position), len(freq)))

    for i, row in enumerate(local_string_position):
        if i < n_samples - samples_per_chunk:
            chunk = local_string_position[i:i+samples_per_chunk]
            abs_fft = 200*np.log10(np.abs(np.fft.fft(chunk, axis=0)))
            output_arr[i,:] = abs_fft
        else:
            output_arr[i,:] = 0


            
    positive_frequency_index = int(len(freq)/2) 
    
    return freq[:positive_frequency_index], output_arr[:,:positive_frequency_index]

def make_animation(A):
    n_frames = np.shape(A)[0]
    n_total_segments = np.shape(A)[1]
    fig = plt.figure()
    maximum = np.max(A)
    minimum = np.min(A)
    padding = 2
    ax = plt.axes(
        xlim=(0 - padding, n_total_segments + padding),
        ylim=(minimum, maximum)
    )

    (line,) = ax.plot([], [], lw=3)
    x = np.arange(0, n_total_segments)

    def plot_A(frame):
        print(frame)
        frame = frame
        line.set_data(x, A[int(frame), :])
        return A[int(frame), :]

    ani = animation.FuncAnimation(fig, plot_A, frames=n_frames, interval=1)
    plt.show()

    return 0;

    

def make_double_animation(A, B, x_b):
    frame_multiplier = 10
    n_total_segments = np.shape(A)[1]
    x1 = np.arange(0, n_total_segments)
    n_frames = np.shape(A)[0]
    n_B_points = np.shape(B)[1]
    maximum_a, minimum_a  = np.max(A), np.min(A)
    maximum_b, minimum_b  = np.max(B), np.min(B)
    maximum_x_b, minimum_x_b  = np.max(x_b), np.min(x_b)
    padding = 2
    fig, (ax1, ax2) = plt.subplots(2,1, sharex=False, sharey=False)
    ax1.set_ylim(bottom=minimum_a, top=maximum_a)
    ax2.set_ylim(bottom=minimum_b, top=maximum_b)
    ax1.set_xlim(left=0, right=n_total_segments)
    ax2.set_xlim(left=minimum_x_b, right=maximum_x_b)

    (line,) = ax1.plot([], [], lw=3)
    (line2,) = ax2.plot([], [], lw=3)
    x_a = np.arange(0, n_total_segments)

    def plot_A(frame):
        frame = frame_multiplier*frame
        frame = frame % n_frames
        line.set_data(x_a, A[int(frame), :])
        line2.set_data(x_b, B[int(frame), :])
        return A[int(frame), :]

    ani = animation.FuncAnimation(fig, plot_A, frames=n_frames, interval=1)
    plt.show()

    return 0;

def bridge_agrafe_index_loc(index_to_phyiscal_position_factor = index_to_phyiscal_position_factor, hammer_location_index = hammer_location_index, side = "bridge"):
    """For comparison with Chaigne 1993-1 p. 6 fig. 3. Converts the location 
    40mm from hammer either to bridge or agrafe side to an actual index on the string.
    side accepted inputs either "bridge" or "agrafe". """
    if side == "bridge":
        shift = -4e-2
    elif side == "agrafe":
        shift = +4e-2
    #print((hammer_location_index * index_to_phyiscal_position_factor + shift) / index_to_phyiscal_position_factor)
    return int((hammer_location_index * index_to_phyiscal_position_factor + shift) / index_to_phyiscal_position_factor), side

def plot_at_string_location(A, string_position_index, side):
    """Plots (x,t) diagram of a point. 
    Noted points are at the bridge and agrafe side both 40mm from the hammer."""
    plt.plot(A[:, string_position_index])
    plt.title(f"String at {side} position, time vs amplitude")
    plt.show()
    return 0

def main():
    plot_index_loc, side = bridge_agrafe_index_loc(side="bridge")
    string_position_arr = simulate_string()
    freqs, fft_arr = make_fft_arr2(string_position_arr, dt, max_t)

    plot_at_string_location(string_position_arr, plot_index_loc, side)
    # make_animation(string_position_arr)
    print(f'{np.shape(string_position_arr)=}')
    print(f'{np.shape(fft_arr)=}')
    make_double_animation(string_position_arr, fft_arr, freqs)


    rate = samp_freq
    data = string_position_arr[:,plot_index_loc]
    scaled = np.int16(data / np.max(np.abs(data)) * 32767)
    write('test.wav', rate, scaled)




if __name__ == "__main__":
    main()
#Implement string results at set points 