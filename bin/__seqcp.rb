#!/usr/bin/env ruby

require 'yaml'
conf = YAML.load_file(ARGV[0])

conf['LIBS'].each do |lib, subconf|
  subconf['FASTQS'].each_index do |i|
    src = subconf['FASTQS'][i]
    cols = src.split(/\/+/)
    dest = "src/#{cols[4]}-#{cols[8]}"
    system "scp -p milou:#{src} #{dest}"
    subconf['FASTQS'][i] = dest
  end
end

outfp = open('conf.converted.yaml', 'w')
YAML.dump(conf, outfp)
outfp.close
