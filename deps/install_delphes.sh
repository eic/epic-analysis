#!/bin/bash
# update or download and install delphes

### download/update source
source environ.sh
if [ -d "$DELPHES_HOME" ]; then
  echo "[+] update Delphes"
  pushd $DELPHES_HOME
  git pull
  popd
else
  echo "[+] downloading Delphes source"
  git clone https://github.com/delphes/delphes.git deps/delphes
fi

### build
echo "[+] building Delphes"
source environ.sh
make delphes-clean
make delphes
