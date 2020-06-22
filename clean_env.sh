cat environment.yml |sed   "s/=\w*$//" | sed "/^prefix/ d" | sed "/^$/ d"  > environment.yml.t
mv environment.yml.t environment.yml
# proper reference to packages installed from git
