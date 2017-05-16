require "bundler/gem_tasks"

task :doc do
  srcs = []
  Dir.glob("ext/numo/linalg/*").each do |d|
    if File.exist?(d+"/extconf.rb")
      sh "cd #{d}; ruby extconf.rb; make clean; make src"
      srcs << d+"/*.c"
    end
  end
  srcs << "lib/numo/linalg.rb"
  sh "yard -m markdown -o yard -r README.md #{srcs.join(' ')}"
end

task :cleandoc do
  sh "rm -r yard .yardoc"
end
