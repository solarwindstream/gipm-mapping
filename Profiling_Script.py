from cProfile import Profile

from pstats import SortKey, Stats



with Profile() as profile:

    import Hann_Window_GIPM

    (Stats(profile).strip_dirs().sort_stats(SortKey.TIME).print_stats())

    # (Stats(profile).strip_dirs().sort_stats(SortKey.CALLS).print_stats())