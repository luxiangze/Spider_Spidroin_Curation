#!/usr/bin/env bash
set -euo pipefail

in=""
out=""

while getopts ":i:o:" opt; do
  case $opt in
    i) in="$OPTARG" ;;
    o) out="$OPTARG" ;;
    *) echo "Usage: $0 -i input.gff3 -o output.gff3" >&2; exit 1 ;;
  esac
done

if [[ -z "${in}" || -z "${out}" ]]; then
  echo "Usage: $0 -i input.gff3 -o output.gff3" >&2
  exit 1
fi

awk -F'\t' -v OFS='\t' '
function get_attr_val(attr,key,  i,a,v){
  split(attr,a,";"); v=""
  for(i in a){ if(a[i] ~ ("^" key "=")){ v=substr(a[i],length(key)+2) } }
  return v
}
function remove_attr(attr,key,  a,i,out){
  split(attr,a,";"); out=""
  for(i=1;i<=length(a);i++){
    if(a[i]!="" && a[i] !~ ("^" key "=")){
      if(out=="") out=a[i]; else out=out";"a[i]
    }
  }
  return out
}
$0 ~ /^#/ { print; next }
$3=="CDS"{
  id = get_attr_val($9,"ID")
  parent = get_attr_val($9,"Parent")
  if(parent==""){ print; next }
  if(id=="") id = parent ":cds"
  key = parent
  c[key]++
  new_id = id ".p" c[key]
  rest = remove_attr($9,"ID")
  if(rest=="") $9="ID=" new_id
  else $9="ID=" new_id ";" rest
  print; next
}
$3=="five_prime_UTR"{
  id = get_attr_val($9,"ID")
  parent = get_attr_val($9,"Parent")
  if(parent==""){ print; next }
  if(id=="") id = parent ":five_prime_UTR"
  key = parent
  c5[key]++
  new_id = id ".p" c5[key]
  rest = remove_attr($9,"ID")
  if(rest=="") $9="ID=" new_id
  else $9="ID=" new_id ";" rest
  print; next
}
$3=="three_prime_UTR"{
  id = get_attr_val($9,"ID")
  parent = get_attr_val($9,"Parent")
  if(parent==""){ print; next }
  if(id=="") id = parent ":three_prime_UTR"
  key = parent
  c3[key]++
  new_id = id ".p" c3[key]
  rest = remove_attr($9,"ID")
  if(rest=="") $9="ID=" new_id
  else $9="ID=" new_id ";" rest
  print; next
}
{ print }
' "$in" > "$out"
