"""
utils.parallel_manager
---------------------
Production-level parallel execution and resource management utilities.

Features:
- ProcessPoolExecutor for CPU-intensive tasks
- ThreadPoolExecutor for I/O-bound tasks
- GPU resource control via CUDA_VISIBLE_DEVICES
- Robust retry with tenacity and exponential backoff
- Structured logging and error handling
- Configurable and reusable for integration
- Real docking and ML task wrappers for Vina, UMol, Boltz2, EquiBind
"""

import os
import logging
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import Callable, List, Any, Optional, Dict, Sequence, Tuple
import tenacity

logger = logging.getLogger("parallel_manager")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s [%(name)s] %(message)s"
)

def robust_retry(
    wait_multiplier: int = 1,
    wait_min: int = 2,
    wait_max: int = 30,
    stop_attempts: int = 5,
    retry_exceptions: Tuple[type, ...] = (Exception,)
) -> Callable:
    """
    Returns a tenacity.retry decorator with exponential backoff.
    """
    return tenacity.retry(
        wait=tenacity.wait_exponential(multiplier=wait_multiplier, min=wait_min, max=wait_max),
        stop=tenacity.stop_after_attempt(stop_attempts),
        retry=tenacity.retry_if_exception_type(retry_exceptions),
        reraise=True
    )

def run_cpu_tasks(
    func: Callable[[Any], Any],
    tasks: Sequence[Any],
    max_workers: Optional[int] = None,
    retry_kwargs: Optional[Dict] = None
) -> List[Any]:
    retry_decorator = robust_retry(**(retry_kwargs or {}))
    wrapped_func = retry_decorator(func)
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_task = {executor.submit(wrapped_func, t): t for t in tasks}
        for future in as_completed(future_to_task):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                logger.error(f"CPU task failed: {e}", exc_info=True)
                raise
    return results

def run_io_tasks(
    func: Callable[[Any], Any],
    tasks: Sequence[Any],
    max_workers: Optional[int] = None,
    retry_kwargs: Optional[Dict] = None
) -> List[Any]:
    retry_decorator = robust_retry(**(retry_kwargs or {}))
    wrapped_func = retry_decorator(func)
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_task = {executor.submit(wrapped_func, t): t for t in tasks}
        for future in as_completed(future_to_task):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                logger.error(f"IO task failed: {e}", exc_info=True)
                raise
    return results

def run_gpu_tasks(
    func: Callable[[Any], Any],
    tasks: Sequence[Any],
    cuda_devices: Sequence[str],
    max_workers: Optional[int] = None,
    retry_kwargs: Optional[Dict] = None
) -> List[Any]:
    retry_decorator = robust_retry(**(retry_kwargs or {}))
    wrapped_func = retry_decorator(func)
    results = []
    max_workers = max_workers or len(cuda_devices)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_task = {
            executor.submit(wrapped_func, t, cuda_device=cuda_devices[i % len(cuda_devices)]): t
            for i, t in enumerate(tasks)
        }
        for future in as_completed(future_to_task):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                logger.error(f"GPU task failed: {e}", exc_info=True)
                raise
    return results

# === REAL TASK FUNCTION WRAPPERS ===
# Import real docking and ML functions from scripts
from scripts.run_docking_multi import run_vina, run_gnina, run_diffdock
from scripts.run_umol import run_umol
from scripts.run_boltz2 import run_boltz2_prediction
from scripts.run_equibind import run_equibind_inference

# Vina wrapper for parallel CPU execution
def vina_task(args):
    """args: (protein, ligand, output_dir, box_params)"""
    return run_vina(*args)

# GNINA wrapper for parallel CPU or GPU execution
def gnina_task(args):
    """args: (protein, ligand, output_dir, use_gpu)"""
    return run_gnina(*args)

# DiffDock wrapper for parallel CPU or GPU execution
def diffdock_task(args):
    """args: (protein, ligand, output_dir)"""
    return run_diffdock(*args)

# UMol wrapper for parallel GPU execution
def umol_task(args, cuda_device=None):
    """args: (protein, ligand, output, ...)"""
    if cuda_device:
        os.environ["CUDA_VISIBLE_DEVICES"] = cuda_device
    return run_umol(*args)

# Boltz2 API prediction wrapper for parallel I/O execution
def boltz2_task(args):
    """args: (protein, ligand, output, ...)"""
    return run_boltz2_prediction(*args)

# EquiBind wrapper for parallel CPU execution
def equibind_task(args):
    """args: (protein, ligand, output_dir, config)"""
    return run_equibind_inference(*args) 