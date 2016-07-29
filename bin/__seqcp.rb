#!/usr/bin/env ruby
require 'yaml'
conf = YAML.load_file(ARGV[0])

system "mkdir -p src/raw"

conf['LIBRARIES'].each do |lib, subconf|
  subconf['FASTQS'].each_index do |i|
    src = subconf['FASTQS'][i]
    cols = src.split(/\/+/)
    dest = "src/raw/#{cols[4]}-#{cols[8]}"
    system "cp --preserve=timestamps #{src} #{dest}"
    subconf['FASTQS'][i] = dest
  end
end

outfp = open('src/conf.converted.yaml', 'w')
YAML.dump(conf, outfp)
outfp.close
