# In a new file: blitzms/parallel/digest_worker.py

import multiprocessing as mp
import time
from multiprocessing import Process, Queue
from queue import Empty
from typing import Any, Literal

from ..proforma.annotation import ProFormaAnnotation
from ..sequence.mod_builder import MOD_BUILDER_INPUT_TYPE, condense_to_peptidoform


class DigestWorker(Process):
    """
    Pre-initialized worker that digests proteins with modifications.

    Avoids serialization overhead by:
    1. Pre-importing all modules in worker process
    2. Processing entire workflow (digest -> modify -> mass -> serialize) in one go
    3. No function pickling - just task parameters
    """

    def __init__(
        self,
        task_queue: Queue,
        result_queue: Queue,
        worker_id: int,
        # Digest parameters (same for all tasks)
        cleave_on: str,
        restrict_before: str,
        restrict_after: str,
        cterminal: bool,
        missed_cleavages: int,
        semi: bool,
        min_len: int | None,
        max_len: int | None,
        # Modification parameters (same for all tasks)
        has_mods: bool,
        nterm_static: MOD_BUILDER_INPUT_TYPE | None,
        cterm_static: MOD_BUILDER_INPUT_TYPE | None,
        internal_static: MOD_BUILDER_INPUT_TYPE | None,
        labile_static: MOD_BUILDER_INPUT_TYPE | None,
        nterm_variable: MOD_BUILDER_INPUT_TYPE | None,
        cterm_variable: MOD_BUILDER_INPUT_TYPE | None,
        internal_variable: MOD_BUILDER_INPUT_TYPE | None,
        labile_variable: MOD_BUILDER_INPUT_TYPE | None,
        max_variable_mods: int,
        use_regex: bool,
        condense_peptidoforms: bool = False,
        invalid_residues: set[str] | None = None,
        use_static_notation: bool = False
    ):
        super().__init__()
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.worker_id = worker_id

        # Store digest parameters
        self.cleave_on = cleave_on
        self.restrict_before = restrict_before
        self.restrict_after = restrict_after
        self.cterminal = cterminal
        self.missed_cleavages = missed_cleavages
        self.semi = semi
        self.min_len = min_len
        self.max_len = max_len

        # Store modification parameters
        self.has_mods = has_mods
        self.nterm_static = nterm_static
        self.cterm_static = cterm_static
        self.internal_static = internal_static
        self.labile_static = labile_static
        self.nterm_variable = nterm_variable
        self.cterm_variable = cterm_variable
        self.internal_variable = internal_variable
        self.labile_variable = labile_variable
        self.max_variable_mods = max_variable_mods
        self.use_regex = use_regex
        self.condense_peptidoforms = condense_peptidoforms
        self.invalid_residues = invalid_residues if invalid_residues is not None else set()
        self.use_static_notation = use_static_notation

    def run(self):
        """Worker main loop - process proteins one at a time."""
        # Import happens ONCE in worker process (not pickled!)
        # These imports are now "warm" for all tasks

        while True:
            try:
                task = self.task_queue.get(timeout=1)

                if task is None:  # Poison pill
                    break

                protein_idx, protein_sequence = task

                try:
                    # Parse protein (happens in worker, no serialization!)
                    protein_annot = ProFormaAnnotation.parse(protein_sequence)

                    # Digest protein
                    peptide_annots: list[ProFormaAnnotation] = list(
                        protein_annot.digest(
                            cleave_on=self.cleave_on,
                            restrict_before=self.restrict_before,
                            restrict_after=self.restrict_after,
                            cterminal=self.cterminal,
                            missed_cleavages=self.missed_cleavages,
                            semi=self.semi,
                            min_len=self.min_len,
                            max_len=self.max_len,
                            return_type="annotation",
                        )
                    )

                    peptide_annots = [p for p in peptide_annots if not any(aa in self.invalid_residues for aa in p.stripped_sequence)]

                    # Apply modifications if needed
                    if self.has_mods:
                        modified_peptide_annots: list[ProFormaAnnotation] = []

                        for peptide in peptide_annots:
                            mod_peptides = list(
                                peptide.build_mods(
                                    nterm_static=self.nterm_static,
                                    cterm_static=self.cterm_static,
                                    internal_static=self.internal_static,
                                    labile_static=self.labile_static,
                                    nterm_variable=self.nterm_variable,
                                    cterm_variable=self.cterm_variable,
                                    internal_variable=self.internal_variable,
                                    labile_variable=self.labile_variable,
                                    max_variable_mods=self.max_variable_mods,
                                    use_regex=self.use_regex,
                                    use_static_notation=self.use_static_notation,
                                    unique_peptidoforms=True,

                                )
                            )
                            modified_peptide_annots.extend(mod_peptides)

                        peptide_annots = modified_peptide_annots


                    serialized_peptide_to_mass: dict[str, float] = {}
                    # ensure that any annotations that have the same mass have different sequences
                    for peptide in peptide_annots:
                        seq = peptide.serialize(include_plus=False, precision=5)
                        if seq not in serialized_peptide_to_mass:
                            serialized_peptide_to_mass[seq] = peptide.mass(precision=5)

                    # Send result back
                    self.result_queue.put(
                        (
                            protein_idx,
                            (list(serialized_peptide_to_mass.keys()), list(serialized_peptide_to_mass.values())),
                            None,  # No error
                        )
                    )

                except Exception as e:
                    # Send error back
                    self.result_queue.put((protein_idx, None, e))

            except Empty:
                continue
            except Exception as e:
                print(f"Worker {self.worker_id} critical error: {e}")
                break


class DigestWorkerPool:
    """Pool of pre-initialized digest workers."""

    def __init__(
        self,
        n_workers: int,
        # Digest parameters
        cleave_on: str,
        restrict_before: str = "",
        restrict_after: str = "",
        cterminal: bool = True,
        missed_cleavages: int = 1,
        semi: bool = False,
        min_len: int | None = 6,
        max_len: int | None = 50,
        # Modification parameters
        has_mods: bool = False,
        nterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
        cterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
        internal_static: MOD_BUILDER_INPUT_TYPE | None = None,
        labile_static: MOD_BUILDER_INPUT_TYPE | None = None,
        nterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
        cterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
        internal_variable: MOD_BUILDER_INPUT_TYPE | None = None,
        labile_variable: MOD_BUILDER_INPUT_TYPE | None = None,
        max_variable_mods: int = 2,
        use_regex: bool = False,
        condense_peptidoforms: bool = False,
        invalid_residues: set[str] | None = None,
        use_static_notation: bool = False,
    ):
        self.n_workers = n_workers
        self.task_queue: Queue = Queue()
        self.result_queue: Queue = Queue()
        self.workers: list[DigestWorker] = []

        # Store parameters for worker initialization
        self.worker_params = {
            "cleave_on": cleave_on,
            "restrict_before": restrict_before,
            "restrict_after": restrict_after,
            "cterminal": cterminal,
            "missed_cleavages": missed_cleavages,
            "semi": semi,
            "min_len": min_len,
            "max_len": max_len,
            "has_mods": has_mods,
            "nterm_static": nterm_static,
            "cterm_static": cterm_static,
            "internal_static": internal_static,
            "labile_static": labile_static,
            "nterm_variable": nterm_variable,
            "cterm_variable": cterm_variable,
            "internal_variable": internal_variable,
            "labile_variable": labile_variable,
            "max_variable_mods": max_variable_mods,
            "use_regex": use_regex,
            "condense_peptidoforms": condense_peptidoforms,
            "invalid_residues": invalid_residues,
            "use_static_notation": use_static_notation,
        }

    def __enter__(self):
        # Start all workers with digest parameters
        for i in range(self.n_workers):
            worker = DigestWorker(
                self.task_queue, self.result_queue, worker_id=i, **self.worker_params
            )
            worker.start()
            self.workers.append(worker)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Send poison pills
        for _ in range(self.n_workers):
            self.task_queue.put(None)

        # Wait for workers
        for worker in self.workers:
            worker.join(timeout=5)
            if worker.is_alive():
                worker.terminate()

        self.workers.clear()

    def digest_proteins(
        self, protein_sequences: list[str]
    ) -> list[tuple[list[str], list[float]]]:
        """
        Digest all proteins and return (peptides, masses) for each.

        Returns:
            List of (serialized_peptides, peptide_masses) tuples, one per protein
        """
        n_proteins = len(protein_sequences)
        
        # Single-threaded mode - bypass multiprocessing overhead
        if self.n_workers == 1:
            start_time = time.time()
            results: list[tuple[list[str], list[float]]] = []
            
            for idx, protein_seq in enumerate(protein_sequences):
                try:
                    peptide_start = time.time()
                    
                    # Parse protein
                    parse_start = time.time()
                    protein_annot = ProFormaAnnotation.parse(protein_seq)
                    parse_time = time.time() - parse_start
                    
                    # Digest protein
                    digest_start = time.time()
                    peptide_annots: list[ProFormaAnnotation] = list(
                        protein_annot.digest(
                            cleave_on=self.worker_params["cleave_on"],
                            restrict_before=self.worker_params["restrict_before"],
                            restrict_after=self.worker_params["restrict_after"],
                            cterminal=self.worker_params["cterminal"],
                            missed_cleavages=self.worker_params["missed_cleavages"],
                            semi=self.worker_params["semi"],
                            min_len=self.worker_params["min_len"],
                            max_len=self.worker_params["max_len"],
                            return_type="annotation",
                        )
                    )
                    digest_time = time.time() - digest_start
                    
                    # Filter invalid residues
                    filter_start = time.time()
                    invalid_residues = self.worker_params["invalid_residues"] or set()
                    peptide_annots = [p for p in peptide_annots if not any(aa in invalid_residues for aa in p.stripped_sequence)]
                    filter_time = time.time() - filter_start
                    
                    # Apply modifications if needed
                    mod_time = 0
                    if self.worker_params["has_mods"]:
                        mod_start = time.time()
                        modified_peptide_annots: list[ProFormaAnnotation] = []
                        
                        for peptide in peptide_annots:
                            mod_peptides = list(
                                peptide.build_mods(
                                    nterm_static=self.worker_params["nterm_static"],
                                    cterm_static=self.worker_params["cterm_static"],
                                    internal_static=self.worker_params["internal_static"],
                                    labile_static=self.worker_params["labile_static"],
                                    nterm_variable=self.worker_params["nterm_variable"],
                                    cterm_variable=self.worker_params["cterm_variable"],
                                    internal_variable=self.worker_params["internal_variable"],
                                    labile_variable=self.worker_params["labile_variable"],
                                    max_variable_mods=self.worker_params["max_variable_mods"],
                                    use_regex=self.worker_params["use_regex"],
                                    use_static_notation=self.worker_params["use_static_notation"],
                                    unique_peptidoforms=True,
                                )
                            )
                            modified_peptide_annots.extend(mod_peptides)
                        
                        peptide_annots = modified_peptide_annots
                        mod_time = time.time() - mod_start
                                        
                    # Serialize and calculate masses
                    serialize_start = time.time()
                    serialized_peptides: list[str] = []
                    peptide_objects: list = []  # Keep references for mass calculation

                    for peptide in peptide_annots:
                        seq = peptide.serialize(include_plus=False, precision=5)
                        if seq not in serialized_peptides:
                            serialized_peptides.append(seq)
                            peptide_objects.append(peptide)

                    serialize_time = time.time() - serialize_start

                    # Calculate masses
                    mass_start = time.time()
                    peptide_masses: list[float] = []
                    for peptide in peptide_objects:
                        peptide_masses.append(peptide.mass(precision=5))
                    mass_time = time.time() - mass_start

                    # Create final dictionary
                    dict_start = time.time()
                    serialized_peptide_to_mass: dict[str, float] = dict(zip(serialized_peptides, peptide_masses))
                    dict_time = time.time() - dict_start

                    peptide_total_time = time.time() - peptide_start

                    # Print timing for first few proteins
                    if idx < 50:
                        print(f"Protein {idx}: parse={parse_time*1000:.1f}ms, digest={digest_time*1000:.1f}ms, "
                            f"filter={filter_time*1000:.1f}ms, mods={mod_time*1000:.1f}ms, "
                            f"serialize={serialize_time*1000:.1f}ms, mass={mass_time*1000:.1f}ms, "
                            f"dict={dict_time*1000:.1f}ms, total={peptide_total_time*1000:.1f}ms, "
                            f"peptides={len(serialized_peptide_to_mass)}")

                    results.append((list(serialized_peptide_to_mass.keys()), list(serialized_peptide_to_mass.values())))
                    
                except Exception as e:
                    raise RuntimeError(f"Error processing protein {idx}: {e}") from e
            
            total_time = time.time() - start_time
            print(f"\nSingle-threaded digest completed: {n_proteins} proteins in {total_time:.2f}s ({n_proteins/total_time:.1f} proteins/sec)")
            return results
        
        # Multi-threaded mode (original code)
        start_time = time.time()
        
        # Submit all proteins
        submit_start = time.time()
        for idx, protein_seq in enumerate(protein_sequences):
            self.task_queue.put((idx, protein_seq))
        submit_time = time.time() - submit_start
        
        # Collect results
        collect_start = time.time()
        results_mt: list[tuple[list[str], list[float]] | None] = [None] * n_proteins
        errors = []
        
        for _ in range(n_proteins):
            protein_idx, result, error = self.result_queue.get()
            
            if error:
                errors.append((protein_idx, error))
            else:
                results_mt[protein_idx] = result
        
        collect_time = time.time() - collect_start
        total_time = time.time() - start_time
        
        print(f"\nMulti-threaded digest completed: {n_proteins} proteins in {total_time:.2f}s")
        print(f"  Submit time: {submit_time:.2f}s, Collect time: {collect_time:.2f}s")
        print(f"  Throughput: {n_proteins/total_time:.1f} proteins/sec")
        
        if errors:
            raise RuntimeError(f"Errors in {len(errors)} proteins: {errors[0][1]}")
        
        return results_mt  # type: ignore
