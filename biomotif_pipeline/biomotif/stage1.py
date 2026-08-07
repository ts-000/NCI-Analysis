from __future__ import annotations

import logging
import math
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .common import min_numeric, safe_float, safe_int, sequence_digest, write_csv, write_json, write_jsonl

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA_URL = "https://data.rcsb.org/rest/v1/core"
FILES_URL = "https://files.rcsb.org/download"
CLUSTER_URL = "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-{identity}.txt"
USER_AGENT = "biomotif-protein-pipeline/1.2"


@dataclass(frozen=True)
class BranchConfig:
    name: str
    resolution_max: float
    rfree_max: float
    min_length: int
    target: int

    def __post_init__(self) -> None:
        if self.name != "protein":
            raise ValueError("Only the protein branch is supported")


ProteinConfig = BranchConfig


@dataclass
class Entity:
    entity_id: str
    is_protein: bool
    sequence: str
    length: int


@dataclass
class Entry:
    branch: str
    pdb_id: str
    title: str
    experimental_methods: list[str]
    resolution_angstrom: float | None
    r_free: float | None
    r_work: float | None
    release_date: str | None
    crystallization_ph: float | None
    crystallization_ph_min: float | None
    crystallization_ph_max: float | None
    crystallization_temperature_K: float | None
    crystallization_details: str
    protein_entity_count: int
    excluded_polymer_entity_count: int
    protein_entity_ids: list[str]
    protein_sequences: list[str]
    protein_lengths: list[int]
    redundancy_key: str
    quality_score: float
    mmcif_path: str | None = None
    structure_factor_path: str | None = None


class RCSBClient:
    def __init__(self, timeout: float = 60.0, retries: int = 5) -> None:
        self.timeout = timeout
        self.retries = retries
        self.local = threading.local()

    def session(self) -> requests.Session:
        session = getattr(self.local, "session", None)
        if session is not None:
            return session

        retry = Retry(
            total=self.retries,
            connect=self.retries,
            read=self.retries,
            status=self.retries,
            backoff_factor=0.7,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset({"GET", "POST"}),
            respect_retry_after_header=True,
            raise_on_status=False,
        )
        adapter = HTTPAdapter(max_retries=retry, pool_connections=20, pool_maxsize=20)
        session = requests.Session()
        session.headers.update({"User-Agent": USER_AGENT})
        session.mount("https://", adapter)
        self.local.session = session
        return session

    def get_json(self, url: str) -> dict[str, Any]:
        response = self.session().get(url, timeout=self.timeout)
        if response.status_code != 200:
            raise RuntimeError(f"GET {url}: HTTP {response.status_code}: {response.text[:500]}")
        data = response.json()
        if not isinstance(data, dict):
            raise RuntimeError(f"GET {url}: expected a JSON object")
        return data

    def download(self, url: str, destination: Path) -> None:
        destination.parent.mkdir(parents=True, exist_ok=True)
        part = destination.with_suffix(destination.suffix + ".part")
        response = self.session().get(url, timeout=self.timeout, stream=True)
        if response.status_code != 200:
            raise RuntimeError(f"Download {url}: HTTP {response.status_code}")
        with part.open("wb") as handle:
            for block in response.iter_content(1024 * 1024):
                if block:
                    handle.write(block)
        if part.stat().st_size == 0:
            part.unlink(missing_ok=True)
            raise RuntimeError(f"Downloaded empty file: {url}")
        part.replace(destination)

    def search(self, cfg: BranchConfig, candidate_pool: int, page_size: int = 500) -> tuple[list[str], dict[str, Any]]:
        query = build_search_payload(cfg, 0, min(page_size, candidate_pool))
        ids: list[str] = []
        start = 0
        while len(ids) < candidate_pool:
            rows = min(page_size, candidate_pool - len(ids))
            query["request_options"]["paginate"] = {"start": start, "rows": rows}
            response = self.session().post(SEARCH_URL, json=query, timeout=self.timeout)
            if response.status_code == 204:
                break
            if response.status_code != 200:
                raise RuntimeError(f"RCSB search HTTP {response.status_code}: {response.text[:1000]}")
            body = response.json()
            result = body.get("result_set", [])
            if not result:
                break
            page = [str(item["identifier"]).upper() for item in result if item.get("identifier")]
            ids.extend(page)
            start += len(page)
            total = body.get("total_count")
            if len(page) < rows or (isinstance(total, int) and start >= total):
                break
        saved = build_search_payload(cfg, 0, candidate_pool)
        return list(dict.fromkeys(ids))[:candidate_pool], saved

    def entry(self, pdb_id: str) -> dict[str, Any]:
        return self.get_json(f"{DATA_URL}/entry/{pdb_id}")

    def polymer_entity(self, pdb_id: str, entity_id: str) -> dict[str, Any]:
        return self.get_json(f"{DATA_URL}/polymer_entity/{pdb_id}/{entity_id}")


def terminal(attribute: str, operator: str, value: Any) -> dict[str, Any]:
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {"attribute": attribute, "operator": operator, "value": value},
    }


def build_search_payload(cfg: BranchConfig, start: int, rows: int) -> dict[str, Any]:
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                terminal("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
                terminal("rcsb_entry_info.resolution_combined", "less_or_equal", cfg.resolution_max),
                terminal("rcsb_entry_info.polymer_entity_count_protein", "greater", 0),
                terminal("rcsb_entry_info.polymer_entity_count_DNA", "equals", 0),
                terminal("rcsb_entry_info.polymer_entity_count_RNA", "equals", 0),
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": start, "rows": rows},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
        },
        "request_info": {"query_id": "biomotif_stage1_protein"},
    }


def classify_entity(data: Mapping[str, Any]) -> Entity:
    poly = data.get("entity_poly") if isinstance(data.get("entity_poly"), Mapping) else {}
    info = data.get("rcsb_polymer_entity") if isinstance(data.get("rcsb_polymer_entity"), Mapping) else {}
    ids = data.get("rcsb_polymer_entity_container_identifiers")
    ids = ids if isinstance(ids, Mapping) else {}
    desc = f"{info.get('rcsb_entity_polymer_type', '')} {poly.get('type', '')}".lower()
    sequence = str(poly.get("pdbx_seq_one_letter_code_can") or poly.get("pdbx_seq_one_letter_code") or "")
    sequence = "".join(sequence.split()).upper()
    return Entity(
        entity_id=str(ids.get("entity_id") or ""),
        is_protein="protein" in desc or "polypeptide" in desc,
        sequence=sequence,
        length=safe_int(poly.get("rcsb_sample_sequence_length"), len(sequence)),
    )


def parse_ph(entry: Mapping[str, Any]) -> tuple[float | None, float | None, float | None, float | None, str]:
    rows = entry.get("exptl_crystal_grow", [])
    if not isinstance(rows, list):
        rows = []
    values: list[float] = []
    lows: list[float] = []
    highs: list[float] = []
    temps: list[float] = []
    details: list[str] = []
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        ph = safe_float(row.get("pH"))
        if ph is not None:
            values.append(ph)
        raw = str(row.get("pdbx_pH_range") or "").replace("–", "-").replace("to", "-")
        nums = [safe_float(tok) for tok in raw.replace(",", " ").replace("-", " ").split()]
        nums = [val for val in nums if val is not None and 0 <= val <= 14]
        if nums:
            lows.append(min(nums))
            highs.append(max(nums))
        temp = safe_float(row.get("temp"))
        if temp is not None:
            temps.append(temp)
        text = str(row.get("pdbx_details") or "").strip()
        if text:
            details.append(text)
    ph = values[0] if values else None
    low = min(lows) if lows else None
    high = max(highs) if highs else None
    if ph is None and low is not None and high is not None:
        ph = 0.5 * (low + high)
    return ph, low, high, min(temps) if temps else None, " | ".join(dict.fromkeys(details))


def parse_entry(cfg: BranchConfig, pdb_id: str, entry: Mapping[str, Any], entities: Sequence[Entity]) -> Entry:
    info = entry.get("rcsb_entry_info") if isinstance(entry.get("rcsb_entry_info"), Mapping) else {}
    struct = entry.get("struct") if isinstance(entry.get("struct"), Mapping) else {}
    accession = entry.get("rcsb_accession_info") if isinstance(entry.get("rcsb_accession_info"), Mapping) else {}
    raw_res = info.get("resolution_combined")
    resolution = min_numeric(raw_res if isinstance(raw_res, list) else [raw_res])
    refine = entry.get("refine", []) if isinstance(entry.get("refine"), list) else []
    r_free = min_numeric(row.get("ls_R_factor_R_free") for row in refine if isinstance(row, Mapping))
    r_work = min_numeric(
        row.get("ls_R_factor_R_work") if row.get("ls_R_factor_R_work") is not None else row.get("ls_R_factor_obs")
        for row in refine if isinstance(row, Mapping)
    )
    methods = [str(row.get("method")) for row in entry.get("exptl", []) if isinstance(row, Mapping) and row.get("method")]
    proteins = [entity for entity in entities if entity.is_protein]
    ph, ph_min, ph_max, temp, details = parse_ph(entry)
    res_score = 0.0 if resolution is None else max(0.0, 1.0 - resolution / cfg.resolution_max)
    r_score = 0.0 if r_free is None else max(0.0, 1.0 - r_free / cfg.rfree_max)
    excluded_count = safe_int(info.get("polymer_entity_count_DNA")) + safe_int(info.get("polymer_entity_count_RNA"))
    return Entry(
        branch="protein",
        pdb_id=pdb_id,
        title=str(struct.get("title") or ""),
        experimental_methods=methods,
        resolution_angstrom=resolution,
        r_free=r_free,
        r_work=r_work,
        release_date=str(accession.get("initial_release_date") or "") or None,
        crystallization_ph=ph,
        crystallization_ph_min=ph_min,
        crystallization_ph_max=ph_max,
        crystallization_temperature_K=temp,
        crystallization_details=details,
        protein_entity_count=safe_int(info.get("polymer_entity_count_protein")),
        excluded_polymer_entity_count=excluded_count,
        protein_entity_ids=[entity.entity_id for entity in proteins],
        protein_sequences=[entity.sequence for entity in proteins],
        protein_lengths=[entity.length for entity in proteins],
        redundancy_key="",
        quality_score=round(100 * (0.65 * res_score + 0.35 * r_score), 3),
    )


def fetch_bundle(client: RCSBClient, cfg: BranchConfig, pdb_id: str) -> Entry:
    entry = client.entry(pdb_id)
    ids = entry.get("rcsb_entry_container_identifiers")
    ids = ids if isinstance(ids, Mapping) else {}
    entity_ids = [str(value) for value in ids.get("polymer_entity_ids", []) or []]
    entities = [classify_entity(client.polymer_entity(pdb_id, entity_id)) for entity_id in entity_ids]
    return parse_entry(cfg, pdb_id, entry, entities)


def passes(entry: Entry, cfg: BranchConfig) -> tuple[bool, str]:
    if "X-RAY DIFFRACTION" not in {method.upper() for method in entry.experimental_methods}:
        return False, "not_xray"
    if entry.resolution_angstrom is None or entry.resolution_angstrom > cfg.resolution_max:
        return False, "resolution"
    if entry.r_free is None:
        return False, "missing_rfree"
    if entry.r_free > cfg.rfree_max:
        return False, "rfree"
    if entry.protein_entity_count < 1 or entry.excluded_polymer_entity_count:
        return False, "composition"
    if not entry.protein_lengths or max(entry.protein_lengths) < cfg.min_length:
        return False, "length"
    return True, "accepted"


def download_cluster_file(client: RCSBClient, identity: int, cache: Path) -> Path:
    path = cache / f"clusters-by-entity-{identity}.txt"
    if not path.exists() or path.stat().st_size == 0:
        client.download(CLUSTER_URL.format(identity=identity), path)
    return path


def cluster_assignments(path: Path, needed: set[str]) -> dict[str, int]:
    found: dict[str, int] = {}
    left = set(needed)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for cluster_id, line in enumerate(handle):
            matches = left.intersection(line.split())
            for key in matches:
                found[key] = cluster_id
            left.difference_update(matches)
            if not left:
                break
    if left:
        logging.warning("%d protein entities were absent from the %d%% cluster file", len(left), int(path.stem.rsplit('-', 1)[-1]))
    return found


def assign_redundancy(entries: Sequence[Entry], clusters: Mapping[str, int]) -> None:
    for entry in entries:
        keys: list[str] = []
        for entity_id in sorted(entry.protein_entity_ids):
            key = f"{entry.pdb_id}_{entity_id}"
            cluster = clusters.get(key)
            keys.append(f"c{cluster}" if cluster is not None else f"exact:{sequence_digest(entry.protein_sequences)}")
        entry.redundancy_key = "protein30:" + ",".join(sorted(keys))


def select_best(entries: Sequence[Entry], target: int) -> tuple[list[Entry], list[Entry]]:
    ordered = sorted(
        entries,
        key=lambda item: (
            item.resolution_angstrom if item.resolution_angstrom is not None else math.inf,
            item.r_free if item.r_free is not None else math.inf,
            -item.quality_score,
            item.pdb_id,
        ),
    )
    picked: list[Entry] = []
    skipped: list[Entry] = []
    seen: set[str] = set()
    for entry in ordered:
        if entry.redundancy_key in seen or len(picked) >= target:
            skipped.append(entry)
            continue
        seen.add(entry.redundancy_key)
        picked.append(entry)
    return picked, skipped


def run_stage1(
    output: Path,
    protein: BranchConfig,
    candidate_pool: int,
    workers: int,
    cluster_identity: int,
    download_structure_factors: bool,
    timeout: float = 60.0,
) -> dict[str, Any]:
    output = output.resolve()
    for name in ("queries", "manifests", "cache", "metadata"):
        (output / name).mkdir(parents=True, exist_ok=True)

    client = RCSBClient(timeout=timeout)
    rejected: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    accepted: list[Entry] = []

    logging.info("Stage 1: searching protein structures")
    ids, query = client.search(protein, candidate_pool)
    write_json(output / "queries" / "protein.json", query)
    logging.info("protein: %d candidates", len(ids))

    with ThreadPoolExecutor(max_workers=workers) as pool:
        jobs = {pool.submit(fetch_bundle, client, protein, pdb_id): pdb_id for pdb_id in ids}
        for done, future in enumerate(as_completed(jobs), 1):
            pdb_id = jobs[future]
            try:
                entry = future.result()
                keep, reason = passes(entry, protein)
                if keep:
                    accepted.append(entry)
                else:
                    rejected.append({"branch": "protein", "pdb_id": pdb_id, "reason": reason})
            except Exception as exc:
                errors.append({"branch": "protein", "pdb_id": pdb_id, "stage": "metadata", "error": str(exc)})
            if done % 100 == 0 or done == len(jobs):
                logging.info("protein metadata: %d/%d", done, len(jobs))

    needed = {f"{entry.pdb_id}_{entity_id}" for entry in accepted for entity_id in entry.protein_entity_ids}
    cluster_file = download_cluster_file(client, cluster_identity, output / "cache")
    assign_redundancy(accepted, cluster_assignments(cluster_file, needed))
    selected, redundant = select_best(accepted, protein.target)

    def download_entry(entry: Entry) -> None:
        cif = output / "protein" / "mmcif" / f"{entry.pdb_id.lower()}.cif.gz"
        if not cif.exists():
            client.download(f"{FILES_URL}/{entry.pdb_id.lower()}.cif.gz", cif)
        entry.mmcif_path = str(cif.relative_to(output))
        if not download_structure_factors:
            return
        sf = output / "protein" / "structure_factors" / f"{entry.pdb_id.lower()}-sf.cif.gz"
        try:
            if not sf.exists():
                client.download(f"{FILES_URL}/{entry.pdb_id.lower()}-sf.cif.gz", sf)
            entry.structure_factor_path = str(sf.relative_to(output))
        except Exception as exc:
            errors.append({"branch": "protein", "pdb_id": entry.pdb_id, "stage": "structure_factors", "error": str(exc)})

    with ThreadPoolExecutor(max_workers=workers) as pool:
        jobs = {pool.submit(download_entry, entry): entry for entry in selected}
        for done, future in enumerate(as_completed(jobs), 1):
            entry = jobs[future]
            try:
                future.result()
            except Exception as exc:
                errors.append({"branch": "protein", "pdb_id": entry.pdb_id, "stage": "coordinates", "error": str(exc)})
            if done % 25 == 0 or done == len(jobs):
                logging.info("coordinate downloads: %d/%d", done, len(jobs))

    downloaded = [entry for entry in selected if entry.mmcif_path]
    summary = {
        "protein": {
            "passed_global_filters": len(accepted),
            "selected_nonredundant": len(selected),
            "target": protein.target,
            "selected_pdb_ids": [entry.pdb_id for entry in selected],
        },
        "protein_selected": len(downloaded),
        "protein_stage2_successful": None,
        "protein_stage2_failed": None,
    }
    rows = [asdict(entry) for entry in downloaded]
    write_csv(output / "manifests" / "selected_structures.csv", rows)
    write_jsonl(output / "manifests" / "selected_structures.jsonl", rows)
    write_csv(output / "manifests" / "redundant_or_unselected.csv", [asdict(entry) for entry in redundant])
    write_csv(output / "manifests" / "rejected.csv", rejected)
    write_json(output / "errors.json", errors)
    write_json(output / "summary.json", summary)
    return {"summary": summary, "errors": errors, "selected": len(downloaded)}
