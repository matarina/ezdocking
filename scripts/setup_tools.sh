#!/usr/bin/env bash
set -euo pipefail

# Helper to fetch P2Rank (prank CLI) and geostd into the repository tree.
# Requires curl, unzip, and git.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
P2RANK_VERSION="2.5.1"
P2RANK_DIR="${ROOT}/p2rank_${P2RANK_VERSION}"
GEOSTD_DIR="${ROOT}/geostd"

echo "Repo root: ${ROOT}"

if [[ ! -d "${P2RANK_DIR}" ]]; then
  echo "Downloading P2Rank ${P2RANK_VERSION}..."
  curl -L "https://github.com/rdk/p2rank/releases/download/${P2RANK_VERSION}/p2rank_${P2RANK_VERSION}.zip" -o "${ROOT}/p2rank_${P2RANK_VERSION}.zip"
  unzip -q "${ROOT}/p2rank_${P2RANK_VERSION}.zip" -d "${ROOT}"
  rm -f "${ROOT}/p2rank_${P2RANK_VERSION}.zip"
else
  echo "P2Rank already present at ${P2RANK_DIR}"
fi

if [[ ! -d "${GEOSTD_DIR}" ]]; then
  echo "Cloning geostd..."
  git clone https://github.com/phenix-project/geostd.git "${GEOSTD_DIR}"
else
  echo "geostd already present at ${GEOSTD_DIR}"
fi

cat <<'EOF'

Setup complete.

Add P2Rank to your PATH for the current shell:
  export PATH="'"${P2RANK_DIR}"'/bin:$PATH"

geostd is cloned to:
  '"${GEOSTD_DIR}"'
The pipeline uses this path by default; no further action needed unless you move it.
EOF
