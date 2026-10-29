# Stream a score file against the ordered reference and write one weight-only
# file per chromosome. The first three score columns must match the reference.
NR == 1 {
  if (NF < 4 || $1 != "SNP" || $2 != "A1" || $3 != "A2") {
    print "Invalid score header" > "/dev/stderr"; failed = 1; exit 1
  }
  header = ""
  for (i = 4; i <= NF; i++) {
    name = $i
    sub(/^SCORE_/, "score_file_" number ".", name)
    header = header (i == 4 ? "" : " ") name
  }
  columns = NF
  next
}
{
  if ((getline refline < map) <= 0) {
    print "Score has more rows than reference" > "/dev/stderr"; failed = 1; exit 1
  }
  split(refline, ref, "\t")
  if (NF != columns || $1 != ref[2] || $2 != ref[3] || $3 != ref[4]) {
    print "Score/reference mismatch at row " NR - 1 > "/dev/stderr"; failed = 1; exit 1
  }
  file = prefix "." ref[1] ".txt"
  if (!(file in opened)) {
    print header > file
    opened[file] = 1
  }
  line = $4
  for (i = 5; i <= NF; i++) line = line " " $i
  print line >> file
}
END {
  if (!failed && NR == 0) {
    print "Empty score file" > "/dev/stderr"; failed = 1
  }
  if (!failed && (getline refline < map) > 0) {
    print "Score has fewer rows than reference" > "/dev/stderr"; failed = 1
  }
  close(map)
  for (file in opened) close(file)
  if (failed) exit 1
}
