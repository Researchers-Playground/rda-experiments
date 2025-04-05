import math
import matplotlib.pyplot as plt

def binary_entropy(epsilon: float) -> float:
    """
    Calculates the binary entropy function h(ε) = -ε log₂(ε) - (1-ε)log₂(1-ε)
    
    Args:
        epsilon: A value between 0 and 1
        
    Returns:
        The binary entropy value
    """
    if epsilon <= 0 or epsilon >= 1:
        raise ValueError("epsilon must be between 0 and 1")
    
    import math
    
    # Handle the case where epsilon is very close to 0 or 1 to avoid log(0)
    if epsilon < 1e-10:
        return -(1-epsilon) * math.log2(1-epsilon)
    if epsilon > 1 - 1e-10:
        return -epsilon * math.log2(epsilon)
    
    return -epsilon * math.log2(epsilon) - (1-epsilon) * math.log2(1-epsilon)

def calculate_rhs(
    delta_sd: float,
    T: int,
    delta_overlap: float,
    delta_overlap_min: float,
    k1: int,
    k2: int,
    epsilon: float,
    N: int
) -> float:
    """
    Calculates the right hand side of the inequality:
    δ ≤ δ_SD + ⌈(T+3)/(Δ_overlap - Δ_overlap,min + 1)⌉ · (k₁2^(h(ε)k₂)e^(-εN/k₁) + k₂e^(-N/k₂))
    
    Args:
        delta_sd: The δ_SD value
        T: Time parameter
        delta_overlap: Δ_overlap value
        delta_overlap_min: Δ_overlap,min value
        k1: k₁ parameter
        k2: k₂ parameter
        epsilon: ε parameter
        N: N parameter
        
    Returns:
        The calculated right hand side value
    """
    import math
    
    # Calculate the ceiling term
    ceiling_term = math.ceil((T + 3) / (delta_overlap - delta_overlap_min + 1))
    
    # Calculate h(ε)
    h_epsilon = binary_entropy(epsilon)
    
    # Calculate the exponential terms
    term1_base2 = pow(2, h_epsilon * k2)
    term1_exp_e = math.exp(-epsilon * N / k1)
    exp_term1 = k1 * term1_base2 * term1_exp_e
    
    exp_term2 = k2 * math.exp(-N / k2)
    
    # Put it all together
    return delta_sd + ceiling_term * (exp_term1 + exp_term2)

def find_max_k1(k2: int, epsilon: float, N: int, target_delta: float = 1e-9) -> int:
    """
    Finds the largest value of k1 that keeps the estimate below target_delta using binary search.
    
    Args:
        k2: k₂ parameter
        epsilon: ε parameter
        N: N parameter
        target_delta: Target upper bound (default: 10^-9)
    
    Returns:
        The largest valid k1 value
    """
    # Time parameters (fixed as per requirements)
    SECONDS_PER_ROUND = 4
    SECONDS_PER_HOUR = 3600
    SECONDS_PER_YEAR = 31557600
    
    T_rounds = int(10 * SECONDS_PER_YEAR / SECONDS_PER_ROUND)
    delta_overlap_rounds = int(6 * SECONDS_PER_HOUR / SECONDS_PER_ROUND)
    delta_overlap_min_rounds = 4
    
    # Binary search for k1
    left, right = 1, 10000  # Reasonable range for k1
    best_k1 = None
    
    while left <= right:
        mid = (left + right) // 2
        result = calculate_rhs(
            delta_sd=0.0,
            T=T_rounds,
            delta_overlap=delta_overlap_rounds,
            delta_overlap_min=delta_overlap_min_rounds,
            k1=mid,
            k2=k2,
            epsilon=epsilon,
            N=N
        )
        
        if result <= target_delta:
            # This k1 works, try to find a larger one
            best_k1 = mid
            left = mid + 1
        else:
            # This k1 is too large
            right = mid - 1
    
    return best_k1

def find_max_k2_with_k1_1(epsilon: float, N: int, target_delta: float = 1e-9) -> int:
    """
    Finds the largest value of k2 that keeps the estimate below target_delta when k1=1 using binary search.
    
    Args:
        epsilon: ε parameter
        N: N parameter
        target_delta: Target upper bound (default: 10^-9)
    
    Returns:
        The largest valid k2 value
    """
    # Time parameters (fixed as per requirements)
    SECONDS_PER_ROUND = 4
    SECONDS_PER_HOUR = 3600
    SECONDS_PER_YEAR = 31557600
    
    T_rounds = int(10 * SECONDS_PER_YEAR / SECONDS_PER_ROUND)
    delta_overlap_rounds = int(6 * SECONDS_PER_HOUR / SECONDS_PER_ROUND)
    delta_overlap_min_rounds = 4
    
    # Binary search for k2
    left, right = 1, 2000  # Increased range for k2
    best_k2 = None
    
    while left <= right:
        mid = (left + right) // 2
        result = calculate_rhs(
            delta_sd=0.0,
            T=T_rounds,
            delta_overlap=delta_overlap_rounds,
            delta_overlap_min=delta_overlap_min_rounds,
            k1=1,  # Fixed to 1
            k2=mid,
            epsilon=epsilon,
            N=N
        )
        
        if result <= target_delta:
            # This k2 works, try to find a larger one
            best_k2 = mid
            left = mid + 1
        else:
            # This k2 is too large
            right = mid - 1
    
    return best_k2

def generate_k1_k2_curve(epsilon: float, N: int, target_delta: float = 1e-9):
    """
    Generates pairs of (k1, k2) values that achieve the target delta.
    For each k2 from max possible down to 1, finds the maximum k1.
    
    Returns:
        List of (k2, k1) pairs
    """
    # First find the maximum k2 possible with k1=1
    max_k2 = find_max_k2_with_k1_1(epsilon, N, target_delta)
    if max_k2 is None:
        return []
    
    curve_points = []
    # For each k2 value from max down to 1
    for k2 in range(max_k2, 0, -1):
        k1 = find_max_k1(k2, epsilon, N, target_delta)
        if k1 is not None:
            curve_points.append((k2, k1))
    
    return curve_points

def calculate_joining_complexity(k1: int, k2: int, n_max: int, t: int = 50) -> float:
    """
    Calculates the total joining complexity using the formula:
    3t + (tn^max)/k₁ + ((t+4)n^max k₂ - 2n^max + (n^max)²)/k₂²
    
    Args:
        k1: k₁ parameter
        k2: k₂ parameter
        n_max: maximum number of nodes
        t: number of bootstrap nodes (default: 50)
    
    Returns:
        The total joining complexity
    """
    term1 = 3 * t
    term2 = (t * n_max) / k1
    term3 = ((t + 4) * n_max * k2 - 2 * n_max + n_max * n_max) / (k2 * k2)
    
    return term1 + term2 + term3

def calculate_get_complexity(k1: int, k2: int, N: int, L_msg: int = 1) -> float:
    """
    Calculates the complexity of the Get operation using the formula:
    (n/(k₁k₂) + 3n²/(k₁k₂²)) · L_msg
    where n = n^max_hon = n^max = N
    
    Args:
        k1: k₁ parameter
        k2: k₂ parameter
        N: maximum number of nodes (n^max)
        L_msg: message length (default: 1)
    
    Returns:
        The Get operation complexity
    """
    n = N  # n = n^max_hon = n^max = N
    term1 = n / (k1 * k2)
    term2 = (3 * n * n) / (k1 * k2 * k2)
    return (term1 + term2) * L_msg

if __name__ == "__main__":
    # Time parameters (for reference)
    SECONDS_PER_ROUND = 4
    SECONDS_PER_HOUR = 3600
    SECONDS_PER_YEAR = 31557600
    
    T_rounds = int(10 * SECONDS_PER_YEAR / SECONDS_PER_ROUND)  # 10 years in rounds
    delta_overlap_rounds = int(6 * SECONDS_PER_HOUR / SECONDS_PER_ROUND)  # 6 hours in rounds
    delta_overlap_min_rounds = 4  # already in rounds
    
    print(f"Time parameters:")
    print(f"T (10 years) = {T_rounds:,} rounds")
    print(f"Δ_overlap (6 hours) = {delta_overlap_rounds:,} rounds")
    print(f"Δ_overlap_min = {delta_overlap_min_rounds} rounds")
    print("Bootstrap nodes (t) = 50")
    print("\n" + "="*50 + "\n")
    
    # Parameters for analysis
    N_values = [1000, 10000, 100000]
    epsilon_values = [0.01, 0.1]
    target_delta = 1e-9
    BOOTSTRAP_NODES = 50  # Fixed number of bootstrap nodes
    L_MSG = 1  # Fixed message length
    
    # Create three subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 8))
    
    # For storing minimum complexity points
    min_join_complexities = []
    min_get_complexities = []
    
    for eps in epsilon_values:
        for N in N_values:
            print(f"\nResults for ε={eps}, N={N}:")
            print("k2\tk1\tJoin Complexity\tGet Complexity")
            print("-" * 50)
            
            # Generate the curve
            curve = generate_k1_k2_curve(eps, N, target_delta)
            
            if curve:  # Only plot if we have points
                k2_values, k1_values = zip(*curve)
                
                # Plot k1-k2 trade-off
                ax1.plot(k2_values, k1_values, marker='o', label=f'ε={eps}, N={N}')
                
                # Calculate and plot join complexities
                join_complexities = [calculate_joining_complexity(k1, k2, N) for k2, k1 in curve]
                ax2.plot(k2_values, join_complexities, marker='o', label=f'ε={eps}, N={N}')
                
                # Calculate and plot get complexities
                get_complexities = [calculate_get_complexity(k1, k2, N) for k2, k1 in curve]
                ax3.plot(k2_values, get_complexities, marker='o', label=f'ε={eps}, N={N}')
                
                # Find minimum complexity points
                min_join_complexity = min(join_complexities)
                min_join_idx = join_complexities.index(min_join_complexity)
                min_join_complexities.append({
                    'epsilon': eps,
                    'N': N,
                    'k1': k1_values[min_join_idx],
                    'k2': k2_values[min_join_idx],
                    'complexity': min_join_complexity
                })
                
                min_get_complexity = min(get_complexities)
                min_get_idx = get_complexities.index(min_get_complexity)
                min_get_complexities.append({
                    'epsilon': eps,
                    'N': N,
                    'k1': k1_values[min_get_idx],
                    'k2': k2_values[min_get_idx],
                    'complexity': min_get_complexity
                })
                
                # Print values
                for k2, k1, join_compl, get_compl in zip(k2_values, k1_values, join_complexities, get_complexities):
                    print(f"{k2}\t{k1}\t{join_compl:.2e}\t{get_compl:.2e}")
            else:
                print("No valid (k2, k1) pairs found for this combination.")
    
    # Configure first subplot (k1-k2 trade-off)
    ax1.set_xlabel('k₂')
    ax1.set_ylabel('k₁')
    ax1.set_title('k₁-k₂ Trade-off Curves')
    ax1.grid(True)
    ax1.legend()
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    
    # Configure second subplot (join complexity)
    ax2.set_xlabel('k₂')
    ax2.set_ylabel('Join Complexity')
    ax2.set_title('Join Complexity vs k₂\n(t=50 bootstrap nodes)')
    ax2.grid(True)
    ax2.legend()
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    
    # Configure third subplot (get complexity)
    ax3.set_xlabel('k₂')
    ax3.set_ylabel('Get Complexity')
    ax3.set_title('Get Complexity vs k₂\n(L_msg=1)')
    ax3.grid(True)
    ax3.legend()
    ax3.set_yscale('log')
    ax3.set_xscale('log')
    
    plt.tight_layout()
    
    # Save the plots
    plt.savefig('complexity_analysis.png')
    plt.close()
    
    # Print minimum complexity points
    print("\nMinimum Join Complexity Points:")
    print("ε\tN\tk1\tk2\tComplexity")
    print("-" * 50)
    for point in min_join_complexities:
        print(f"{point['epsilon']}\t{point['N']}\t{point['k1']}\t{point['k2']}\t{point['complexity']:.2e}")
        
    print("\nMinimum Get Complexity Points:")
    print("ε\tN\tk1\tk2\tComplexity")
    print("-" * 50)
    for point in min_get_complexities:
        print(f"{point['epsilon']}\t{point['N']}\t{point['k1']}\t{point['k2']}\t{point['complexity']:.2e}")