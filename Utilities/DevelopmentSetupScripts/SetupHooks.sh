#!/usr/bin/env bash
#==========================================================================
#
#   Copyright Insight Software Consortium
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#          http://www.apache.org/licenses/LICENSE-2.0.txt
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
#==========================================================================*/


# Run this script to set up the git hooks for committing changes to ANTs.

egrep-q() {
  egrep "$@" >/dev/null 2>/dev/null
}

die() {
  echo 'failure during hook setup' 1>&2
  echo '-------------------------' 1>&2
  echo '' 1>&2
  echo "$@" 1>&2
  exit 1
}

u=$(cd "$(echo "$0"|sed 's/[^/]*$//')"; pwd)
cd "$u/../../.git/hooks"

# Symlink the hooks in from a local dir.
for file in ../../Utilities/Hooks/* ../../Utilities/Hooks/.*; do
  base=$(basename "$file")
  if [ x"$base" == x"." -o x"$base" == x".." ]; then
    continue
  fi
  rm "$base"
  ln -s "$file"
done

# Set up uncrustify hook.
echo "Setting up the uncrustify hook..."
git config hooks.uncrustify.conf "Utilities/Maintenance/uncrustify.cfg"

# Set up KWStyle hook.
echo "Setting up the KWStyle hook..."
git config hooks.KWStyle.conf "Utilities/KWStyle/kws.xml.in"
git config hooks.KWStyle.overwriteRulesConf "Utilities/KWStyle/Overwrite.txt"
git config hooks.KWStyle false

# Set up cppcheck hook.
echo "Setting up the cppcheck hook..."
git config hooks.cppcheck true

# Set up hook chaining.
echo "Setting up hook chaining: prepare-commit-msg, commit-msg, pre-commit"
git config hooks.chain-prepare-commit-msg Utilities/Hooks/prepare-commit-msg
git config hooks.chain-commit-msg Utilities/Hooks/commit-msg
git config hooks.chain-pre-commit Utilities/Hooks/pre-commit

echo "Done."
