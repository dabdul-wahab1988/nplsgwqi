from __future__ import annotations

from reviewer_revision_artifacts import (
    generate_irrigation_disagreement,
    generate_risk_uncertainty_decomposition,
    generate_wqi_benchmarking,
)


def main() -> None:
    print("Generating WQI benchmarking + weight uncertainty…")
    generate_wqi_benchmarking()

    print("Generating irrigation disagreement tables…")
    generate_irrigation_disagreement()

    print("Generating risk uncertainty decomposition…")
    generate_risk_uncertainty_decomposition()

    print("OK: stage2 artifacts generated.")


if __name__ == "__main__":
    main()

