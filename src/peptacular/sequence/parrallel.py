import multiprocessing as mp
import sys
from functools import partial
from multiprocessing.pool import Pool, ThreadPool
from typing import Any, Literal, TypeVar
from collections.abc import Callable, Sequence
import atexit

from ..constants import ParrallelMethod, ParrallelMethodLiteral


T = TypeVar("T")

# Global pool cache
_pool_cache: dict[tuple[str, int], Pool | ThreadPool] = {}
_pool_lock = mp.Lock()


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
    method: ParrallelMethod, n_workers: int, reuse_pool: bool = True
) -> Pool | ThreadPool:
    """
    Get an existing pool or create a new one.

    :param method: Process or thread pool
    :param n_workers: Number of workers
    :param reuse_pool: If True, reuse existing pools
    :return: Pool instance
    """
    if not reuse_pool:
        if method == ParrallelMethod.THREAD:
            return ThreadPool(processes=n_workers)
        else:
            return Pool(processes=n_workers)

    pool_key = (method.value, n_workers)

    with _pool_lock:
        if pool_key not in _pool_cache:
            if method == ParrallelMethod.THREAD:
                _pool_cache[pool_key] = ThreadPool(processes=n_workers)
            else:
                _pool_cache[pool_key] = Pool(processes=n_workers)

        return _pool_cache[pool_key]


def cleanup_pools():
    """
    Clean up all cached pools. Call this at program exit if needed.
    """
    global _pool_cache
    with _pool_lock:
        for pool in _pool_cache.values():
            try:
                pool.close()
                pool.join()
            except Exception:
                pass
        _pool_cache.clear()


def parallel_apply_internal(
    func: Callable[..., T],
    items: Sequence[Any],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
    reuse_pool: bool = True,
    verbose: bool = False,
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
    :param reuse_pool: If True (default), reuse existing pools for better performance
    :param verbose: If True, print worker and chunksize info
    :param func_kwargs: Keyword arguments to pass to the function
    :return: List of results in the same order as input items
    """
    # Convert to list if needed
    items_list = list(items)

    # Handle empty input
    if not items_list:
        return []

    method_enum = (
        ParrallelMethod(method) if method is not None else ParrallelMethod.PROCESS
    )

    # Handle sequential execution
    if method_enum == ParrallelMethod.SEQUENTIAL:
        return [func(item, **func_kwargs) for item in items_list]

    # Handle parallel execution
    if n_workers is None:
        n_workers = mp.cpu_count()

    if verbose:
        print(f"Using n_workers: {n_workers}")

    # Optimize chunksize calculation
    if chunksize is None:
        # Better heuristic: balance overhead vs parallelism
        total_items = len(items_list)
        if total_items < n_workers:
            chunksize = 1
        else:
            # Aim for ~4 chunks per worker to balance load
            chunksize = max(1, total_items // (n_workers * 4))

    if verbose:
        print(f"Using chunksize: {chunksize}")

    wrapper = partial(_apply_wrapper, func=func, func_kwargs=func_kwargs)

    # Get or create pool
    pool = _get_or_create_pool(method_enum, n_workers, reuse_pool)

    if reuse_pool:
        # Don't use context manager - pool is managed globally
        return pool.map(wrapper, items_list, chunksize=chunksize)
    else:
        # Create temporary pool
        with pool:
            return pool.map(wrapper, items_list, chunksize=chunksize)


atexit.register(cleanup_pools)
