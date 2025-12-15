import multiprocessing as mp
import sys
from functools import partial
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool
from typing import Any, Callable, Literal, Sequence, TypeVar

from ..constants import ParrallelMethod, ParrallelMethodLiteral

T = TypeVar("T")


def _is_gil_disabled() -> bool:
    """
    Detect if the GIL is disabled (Python 3.13+ with free-threading).
    :return: True if GIL is disabled, False otherwise
    """
    # Python 3.13+ with free-threading
    if sys.version_info >= (3, 13):
        try:
            if hasattr(sys, "_is_gil_enabled"):
                result = not sys._is_gil_enabled()
                return result
        except AttributeError:
            pass
    return False


def _get_optimal_method() -> Literal["process", "thread"]:
    """
    Automatically determine the best parallelization method.

    :return: 'thread' if GIL is disabled, 'process' otherwise
    """
    return "thread" if _is_gil_disabled() else "process"


def _apply_wrapper(item: Any, func: Callable[..., T], func_kwargs: dict[str, Any]) -> T:
    """
    Wrapper function at module level so it can be pickled.

    :param item: Single item to process
    :param func: Function to apply
    :param func_kwargs: Keyword arguments for the function
    :return: Result of func(item, **func_kwargs)
    """
    return func(item, **func_kwargs)


def parallel_apply_internal(
    func: Callable[..., T],
    items: Sequence[Any],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
    **func_kwargs: Any,
) -> list[T]:
    """
    Internal function for parallel processing.
    Automatically detects if GIL is disabled and uses threading for better performance.
    Otherwise defaults to multiprocessing.
    Results are returned in the same order as the input items.

    :param func: Function to apply to each item
    :param items: Sequence of items to process
    :param n_workers: Number of worker processes/threads. If None, uses CPU count
    :param chunksize: Number of items per chunk. If None, auto-calculated
    :param method: 'process', 'thread', 'sequential', or None (auto-detect). Default is None.
    :param func_kwargs: Keyword arguments to pass to the function
    :return: List of results in the same order as input items
    """
    # Convert to list if needed
    items_list = list(items)

    method_enum = (
        ParrallelMethod(method) if method is not None else ParrallelMethod.PROCESS
    )

    # Handle sequential execution
    if method_enum == ParrallelMethod.SEQUENTIAL:
        return [func(item, **func_kwargs) for item in items_list]

    # Handle parallel execution
    if n_workers is None:
        n_workers = mp.cpu_count()

    if chunksize is None:
        chunksize = max(1, len(items_list) // (n_workers * 2))
        print(f"Using chunksize: {chunksize}")

    wrapper = partial(_apply_wrapper, func=func, func_kwargs=func_kwargs)

    # Create a new pool for each call
    PoolClass = ThreadPool if method_enum == ParrallelMethod.THREAD else Pool
    with PoolClass(processes=n_workers) as pool:
        return pool.map(wrapper, items_list, chunksize=chunksize)
