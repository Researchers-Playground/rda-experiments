from simulation_core import *
import matplotlib.pyplot as plt

if __name__ == "__main__":
    cols = 100
    n_total = 2500  # say 10000 nodes and 2500 are honest
    n_init = 20
    n_warmup = n_total - n_init
    steps = 5000
    churn = 50
    lifetime_per_party = n_total // churn

    # Range of row values to iterate over
    row_range = [1, 5, 10, 25]

    # Define different markers and colors for each plot
    markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'h']
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple', 'brown']

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))  # Two subplots side by side
    fig.suptitle(f'Protocol Simulation with {cols} columns, {n_total} honest parties, staying for {lifetime_per_party} steps', fontsize=14)

    for i, rows in enumerate(row_range):
        params = ProtocolParameters(k1=rows, k2=cols, delta_sub=1, m=100)
        schedule = generate_schedule(n_init=n_init, n_warmup=n_warmup, churn=churn, steps=steps)
        statistics = simulate_protocol_run(schedule, params)

        # Left plot - Corruption Graph
        # Make relative first
        corruption_graph_relative = [100 * x / cols for x in statistics.corruption_graph]
        axes[0].plot(corruption_graph_relative, label=f'Rows = {rows}',
                     marker=markers[i % len(markers)], color=colors[i % len(colors)], markevery=500)

        # Right plot - Connection Graph
        axes[1].plot(statistics.connections_graph, label=f'Rows = {rows}',
                     marker=markers[i % len(markers)], color=colors[i % len(colors)], markevery=500)

    # Customize left plot
    axes[0].set_xlabel("Time Steps")
    axes[0].set_ylabel("Max Fraction of Corrupted Symbols (in %)")
    axes[0].set_title('Corruption Graphs')
    axes[0].set_ylim(0, 100)
    axes[0].legend()
    axes[0].grid(True)

    # Customize right plot
    axes[1].set_xlabel("Time Steps")
    axes[1].set_ylabel("Max Number of Peers")
    axes[1].set_title('Connection Graphs')
    axes[1].legend()
    axes[1].grid(True)

    plt.tight_layout()
    plt.show()
