##
## Step 3d - merge the counts
##

step3d_sources = Array.new
LIBWELLIDS.each do |libwellid|
  step3d_sources.push("tmp/cg/#{libwellid}.step3b")
end

file 'out/cg/reads_all.txt.gz' => step3d_sources do |t|
  outfp = open("| gzip -c > #{t.name}", 'w')
  outfp.puts (['Gene'] + LIBWELLIDS).join("\t")

  tmp = Array.new
  LIBWELLIDS.each do |libwellid|
    tmp.push("tmp/cg/#{libwellid}.step3b")
  end
  infp = open("| cat "+tmp.join(" | join -t '\t' -j 1 - "))
  while line = infp.gets
    outfp.puts line
  end
  infp.close
  outfp.close
end

task :clean_step3d do
  sh 'rm out/cg/reads_all.txt.gz'
end
