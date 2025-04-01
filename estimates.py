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
    print("\n" + "="*50 + "\n")
    
    # Parameters for k1-k2 analysis
    N_values = [1000, 10000, 100000]
    epsilon_values = [0.01, 0.1]
    target_delta = 1e-9
    
    # Generate and display tables and plots for each combination
    plt.figure(figsize=(15, 10))
    
    for eps in epsilon_values:
        for N in N_values:
            print(f"\nResults for ε={eps}, N={N}:")
            print("k2\tk1")
            print("-" * 20)
            
            # Generate the curve
            curve = generate_k1_k2_curve(eps, N, target_delta)
            
            # Print the table
            #for k2, k1 in curve:
            #    print(f"{k2}\t{k1}")
            
            # Plot the curve
            if curve:  # Only plot if we have points
                k2_values, k1_values = zip(*curve)
                plt.plot(k2_values, k1_values, marker='o', label=f'ε={eps}, N={N}')
            else:
                print("No valid (k2, k1) pairs found for this combination.")    
    plt.xlabel('k₂')
    plt.ylabel('k₁')
    plt.title('k₁-k₂ Trade-off Curves for Different ε and N Values\n(δ ≤ 10⁻⁹)')
    plt.grid(True)
    plt.legend()
    plt.yscale('log')
    plt.xscale('log')
    
    # Save the plot
    plt.savefig('k1_k2_tradeoff.png')
    plt.close()