import atexit
import multiprocessing as mp
import sys
from functools import partial
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool
from typing import Any, Callable, Literal, Sequence, TypeVar, Union

from ..constants import ParrallelMethod, ParrallelMethodLiteral

T = TypeVar("T")

# Global pool cache: (method, n_workers) -> pool instance
_POOL_CACHE = {}


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


def _get_or_create_pool(
    method: Literal["process", "thread"], n_workers: int
) -> Union[Pool, ThreadPool]:
    """
    Get a cached pool or create a new one if it doesn't exist.

    :param method: 'process' or 'thread'
    :param n_workers: Number of workers
    :return: Pool or ThreadPool instance
    """
    cache_key = (method, n_workers)

    if cache_key not in _POOL_CACHE:
        PoolClass = ThreadPool if method == "thread" else Pool
        _POOL_CACHE[cache_key] = PoolClass(processes=n_workers)

    return _POOL_CACHE[cache_key]


def _cleanup_pools() -> None:
    """
    Clean up all cached pools. Called automatically at exit.
    """
    for pool in _POOL_CACHE.values():
        pool.close()
        pool.join()
    _POOL_CACHE.clear()


def clear_pool_cache() -> None:
    """
    Manually clear the pool cache and terminate all workers.
    Useful for explicit cleanup or when you want to force recreation of pools.
    """
    _cleanup_pools()


def get_pool_cache_info() -> dict[str, Any]:
    """
    Get information about the current pool cache.
    
    :return: Dictionary with cache statistics:
        - 'size': Number of cached pools
        - 'pools': List of (method, n_workers) tuples for each cached pool
    """
    return {
        "size": len(_POOL_CACHE),
        "pools": list(_POOL_CACHE.keys()),
    }


# Register cleanup function to run at exit
atexit.register(_cleanup_pools)


def parallel_apply_internal(
    func: Callable[..., T],
    items: Sequence[Any],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
    use_cache: bool = True,
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
    :param use_cache: If True, reuse cached pools. If False, create a new pool each time. Default is True.
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

    # Use cached pool or create new one
    if use_cache:
        pool_method = "thread" if method_enum == ParrallelMethod.THREAD else "process"
        pool = _get_or_create_pool(pool_method, n_workers)
        return pool.map(wrapper, items_list, chunksize=chunksize)
    else:
        # Create a new pool for this call only (old behavior)
        PoolClass = ThreadPool if method_enum == ParrallelMethod.THREAD else Pool
        with PoolClass(processes=n_workers) as pool:
            return pool.map(wrapper, items_list, chunksize=chunksize)
