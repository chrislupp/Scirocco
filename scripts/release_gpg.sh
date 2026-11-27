#!/usr/bin/env bash
# Generate a one-day GPG key and configure git to sign commits and tags as
# github-actions[bot]. The release workflow runs this in each job that
# writes commits or tags; the key exists only on that runner.
set -euo pipefail
params="$(mktemp)"
cat > "$params" <<GPGEOF
%no-protection
Key-Type: RSA
Key-Length: 4096
Subkey-Type: RSA
Subkey-Length: 4096
Name-Real: github-actions[bot]
Name-Email: 41898282+github-actions[bot]@users.noreply.github.com
Expire-Date: 1d
GPGEOF
gpg --batch --gen-key "$params"
rm -f "$params"
key_id=$(gpg --list-secret-keys --keyid-format=long | grep sec | awk '{print $2}' | cut -d'/' -f2)

git config user.name  "github-actions[bot]"
git config user.email "41898282+github-actions[bot]@users.noreply.github.com"
git config user.signingkey "${key_id}"
git config commit.gpgsign true
git config tag.gpgsign true
