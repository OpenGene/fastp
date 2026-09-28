#!/usr/bin/env bash
# Create or update the benchmark comment on a PR.
# usage: sticky_comment.sh OWNER/REPO PR_NUMBER BODY_FILE   (needs GH_TOKEN with pull-requests: write)
set -euo pipefail
repo=$1 pr=$2 body=$3
marker='<!-- fastp-benchmark -->'
[[ $pr =~ ^[0-9]+$ ]] || { echo "bad PR number: $pr" >&2; exit 1; }
[[ $(head -c ${#marker} "$body") == "$marker" ]] || { echo "body lacks the benchmark marker" >&2; exit 1; }
head -c 60000 "$body" > body.trimmed  # GitHub's comment limit is 65536
id=$(gh api "repos/$repo/issues/$pr/comments" --paginate \
  --jq ".[] | select(.user.login == \"github-actions[bot]\" and (.body | startswith(\"$marker\"))) | .id" | tail -1)
if [[ -n $id ]]; then
  gh api -X PATCH "repos/$repo/issues/comments/$id" -F body=@body.trimmed >/dev/null
else
  gh api "repos/$repo/issues/$pr/comments" -F body=@body.trimmed >/dev/null
fi
