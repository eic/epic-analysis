#!/usr/bin/env ruby

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

# truncate a config file to the specified number of ROOT files

if ARGV.length < 2
  $stderr.puts "USAGE: #{$0} [config_file] [num_files]"
  exit 2
end
configN = ARGV[0]
num     = ARGV[1].to_i

outN = configN + ".truncated"
out = File.open outN, 'w'
cnt = 0

File.readlines(configN).each do |line|
  if line.match? /\.root/
    unless line.match? /^[#:]/
      cnt += 1
    end
  end
  out.puts line
  break if cnt >= num
end

out.close
system "cat #{outN}"
puts "#{'-'*40}\nwrote #{outN}"
