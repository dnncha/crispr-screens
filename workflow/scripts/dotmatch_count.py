import csv
import json
import subprocess
import tempfile
from pathlib import Path


"""Run DotMatch and adapt its per-target table to this workflow's count input."""

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with tempfile.TemporaryDirectory(prefix="crispr-screens-dotmatch-") as temporary:
    raw_counts = Path(temporary) / "counts.tsv"
    command = [
        "dotmatch",
        "count",
        "--targets",
        str(snakemake.input.targets),
        "--reads",
        str(snakemake.input.fq),
        "--target-start",
        str(snakemake.params.target_start),
        "--target-length",
        str(snakemake.params.target_length),
        "--k",
        str(snakemake.params.k),
        "--metric",
        "hamming",
        "--ambiguity-policy",
        str(snakemake.params.ambiguity_policy),
        "--out",
        str(raw_counts),
        "--summary",
        str(snakemake.output.summary),
    ]
    completed = subprocess.run(command, capture_output=True, text=True)

    with log_path.open("w", encoding="utf-8") as log:
        if completed.stdout:
            log.write(completed.stdout)
            if not completed.stdout.endswith("\n"):
                log.write("\n")
        if completed.stderr:
            log.write(completed.stderr)
            if not completed.stderr.endswith("\n"):
                log.write("\n")
        if completed.returncode != 0:
            raise subprocess.CalledProcessError(
                completed.returncode,
                command,
                output=completed.stdout,
                stderr=completed.stderr,
            )

        summary = json.loads(Path(snakemake.output.summary).read_text(encoding="utf-8"))
        samples = summary.get("samples", [])
        if len(samples) > 1:
            raise ValueError("DotMatch count summary must contain exactly one sample")
        sample_summary = samples[0] if samples else summary
        total_reads = int(sample_summary["total_reads"])
        assigned_unique = int(sample_summary["assigned_unique"])
        assignment_rate = 100 * assigned_unique / total_reads if total_reads else 0
        log.write(f"DotMatch unique assignment rate: {assignment_rate:.2f}%\n")
        log.write(f"DotMatch total reads: {total_reads}\n")
        log.write(f"DotMatch unique assignments: {assigned_unique}\n")
        log.write(f"DotMatch ambiguous reads: {int(sample_summary['ambiguous'])}\n")
        log.write(f"DotMatch unmatched reads: {int(sample_summary['unmatched'])}\n")
        log.write(f"DotMatch invalid reads: {int(sample_summary['invalid'])}\n")

    with raw_counts.open(encoding="utf-8", newline="") as source, Path(
        snakemake.output.counts
    ).open("w", encoding="utf-8") as target:
        for row in csv.DictReader(source, delimiter="\t"):
            count_column = next(
                (name for name in row if name.endswith("_count_total")),
                "count_total",
            )
            if count_column not in row:
                raise ValueError("DotMatch count output has no total-count column")
            target.write(f"{row[count_column]} {row['target_id']}\n")
