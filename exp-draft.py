from typing import List, NamedTuple, Dict, Tuple, Optional
import random
import matplotlib.pyplot as plt

class PartyGrid:
    def __init__(self, k1: int, k2: int, seed: Optional[int] = None):
        """A k1 x k2 grid where each cell contains a list of integers (= parties)."""
        self.k1 = k1
        self.k2 = k2
        self.grid = [[[] for _ in range(k2)] for _ in range(k1)]  # Grid of lists
        self.party_to_cell: Dict[int, Tuple[int, int]] = {}  # Maps party to cell coordinates
        self.rng = random.Random(seed)  # Initialize random generator with seed

    def add_party(self, party: int) -> bool:
        """Assigns a party to a random cell. Returns True if successful, False if the party already exists."""
        if party in self.party_to_cell:
            return False  # Party already exists in a cell
        row, col = self.rng.randint(0, self.k1 - 1), self.rng.randint(0, self.k2 - 1)
        self.grid[row][col].append(party)
        self.party_to_cell[party] = (row, col)
        return True

    def get_cell_of_party(self, party: int) -> Optional[Tuple[int, int]]:
        """Returns the coordinates of the cell containing the given party, or None if not found."""
        return self.party_to_cell.get(party)

    def get_parties_in_cell(self, row: int, col: int) -> List[int]:
        """Returns a list of parties in the specified cell."""
        if 0 <= row < self.k1 and 0 <= col < self.k2:
            return self.grid[row][col]
        return []  # Invalid cell coordinates

    def get_parties_in_row(self, row: int) -> List[int]:
        """Returns all parties in the specified row."""
        if 0 <= row < self.k1:
            return [party for cell in self.grid[row] for party in cell]
        return []  # Invalid row index

    def get_parties_in_column(self, col: int) -> List[int]:
        """Returns all parties in the specified column."""
        if 0 <= col < self.k2:
            return [party for row in self.grid for party in row[col]]
        return []  # Invalid column index

    def remove_party(self, party: int) -> bool:
        """Removes a party from its cell. Returns True if successful, False if the party was not found."""
        if party in self.party_to_cell:
            row, col = self.party_to_cell.pop(party)
            self.grid[row][col].remove(party)
            return True
        return False

    def get_corrupted_columns_for_row(self, row: int) -> int:
        """
        Returns the number of columns that are corrupted for a row.
        """

        corrupted_cols = [c for c in range(self.k2) if len(self.grid[row][c]) == 0]
        return len(corrupted_cols)

    def get_max_corrupted_columns(self) -> int:
        """
        Returns max_P # corrupted columns for P
        """
        return max((self.get_corrupted_columns_for_row(row) for row in range(self.k1)), default=0)

class ProtocolParameters(NamedTuple):
    k1: int
    k2: int
    delta_sub: int
    m: int

class JoinEvent(NamedTuple):
    party: int

class LeaveEvent(NamedTuple):
    party: int

Event = NamedTuple("Event", [
    ("type", str),
    ("data", JoinEvent or LeaveEvent)
])

def generate_schedule(n_init: int, D: int, steps: int) -> List[List[Event]]:
    """Generates a join-leave schedule."""
    schedule = []
    party_counter = 0
    active_parties = []

    # Initial join events
    init_events = [Event("join", JoinEvent(party_counter + i)) for i in range(n_init)]
    party_counter += n_init
    active_parties.extend(range(n_init))
    schedule.append(init_events)

    # Subsequent time steps
    for _ in range(steps):
        step_events = []
        # D parties join
        for _ in range(D):
            step_events.append(Event("join", JoinEvent(party_counter)))
            active_parties.append(party_counter)
            party_counter += 1
        # D parties leave (FIFO rule)
        for _ in range(D):
            if active_parties:
                step_events.append(Event("leave", LeaveEvent(active_parties.pop(0))))
        schedule.append(step_events)

    return schedule

def simulate_protocol_run(schedule: List[List[Event]], protocol_params: ProtocolParameters) -> List[int]:
    """
    Simulates the execution of the protocol based on the given schedule and protocol parameters.
    """
    assert len(schedule) > 0, "schedule length should be non-zero"

    # data structure we need
    grid = PartyGrid(protocol_params.k1, protocol_params.k2, seed=42)

    # initialization (time tau = 0)
    init_events = schedule[0]
    assert all(event.type == "join" for event in init_events), "All initial events must be join events"
    for event in init_events:
        grid.add_party(event.data.party)
    corruption_graph = [grid.get_max_corrupted_columns()]

    # times tau > 0, protocol run
    for events in schedule[1:]:
        # execute join and leave events
        for event in events:
            if event.type == "join":
                grid.add_party(event.data.party)
            elif event.type == "leave":
                grid.remove_party(event.data.party)
        corruption_graph.append(grid.get_max_corrupted_columns())

    return corruption_graph

if __name__ == "__main__":
    rows = 10
    cols = 10
    params = ProtocolParameters(k1=rows, k2=cols, delta_sub=1, m=100)
    schedule = generate_schedule(n_init=200, D=50, steps=20)
    corruption_graph = simulate_protocol_run(schedule, params)

    # Plot the graph
    plt.plot(corruption_graph)
    plt.xlabel("Time Steps")
    plt.ylabel("Max Corrupted Columns")
    plt.title("Protocol Run Simulation")
    plt.ylim(0, cols)

    # Show the plot
    plt.show()
