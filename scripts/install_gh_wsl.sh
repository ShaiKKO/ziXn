#!/usr/bin/env bash
set -euo pipefail

# Requires root. Use: sudo -n bash scripts/install_gh_wsl.sh

if ! command -v apt-get >/dev/null 2>&1; then
  echo "apt-get not found. Unsupported WSL distro." >&2
  exit 1
fi

export DEBIAN_FRONTEND=noninteractive
apt-get update -y
apt-get install -y curl ca-certificates gnupg lsb-release

# Install GitHub CLI official repo key
curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg \
  | dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg status=none
chmod go+r /usr/share/keyrings/githubcli-archive-keyring.gpg

arch="$(dpkg --print-architecture)"
dist="stable"
echo "deb [arch=${arch} signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages ${dist} main" \
  > /etc/apt/sources.list.d/github-cli.list

apt-get update -y
apt-get install -y gh

echo "Installed:" $(gh --version | head -n1)


