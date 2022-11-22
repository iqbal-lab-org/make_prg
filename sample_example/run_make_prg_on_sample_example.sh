#!/usr/bin/env bash
set -eu

# configs
version="0.4.0"
make_prg="./make_prg_${version}"
# TODO: update to release URL when we have it
make_prg_URL="https://www.dropbox.com/s/jfzuxjkjd5ke15x/make_prg_0.4.0?dl=1"


if [ ! -f ${make_prg} ]; then
  echo "${make_prg} not found, downloading it..."
  wget ${make_prg_URL} -O ${make_prg}
  chmod +x ${make_prg}
fi

echo "Building PRGs from MSAs..."
from_msa_command_line="${make_prg} from_msa --input msas/ --output-prefix msas_output/sample --force"
echo "Running ${from_msa_command_line}"
${from_msa_command_line}

echo "Updating PRGs with denovo paths..."
update_command_line="${make_prg} update --update-DS msas_output/sample.update_DS.zip --denovo-paths denovo_paths/denovo_paths.txt --output-prefix msas_updated/updated_sample --force"
echo "Running ${update_command_line}"
${update_command_line}

echo "Sample example ran successfully!"
