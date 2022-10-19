#!/usr/bin/env ruby
# remove any partial, final event from a HEPMC file
# (not efficient, this is only used for CI tests)

if ARGV.length<2
  $stderr.puts "USAGE: #{$0} [hepmc file] [output file]"
  exit 2
end
streamFile = ARGV[0]
outFile    = ARGV[1]

puts "TRIMMING #{streamFile}"

out = File.open outFile, 'w'
buffer = []
File.readlines(streamFile).each do |line|
  if line.match? /^E /
    buffer.each do |output| out.puts output end
    buffer.clear
  end
  buffer << line
end
out.close

puts " -> #{outFile}"
