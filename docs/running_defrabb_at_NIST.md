# Running DeFrABB using NIST Compute Resources (NIST Internal Use Only)
The following usage documentation is specific for running on NIST systems. 
Please use these as a starting point for running the pipeline locally. Contact Nathan Olson at nolson@nist.gov, with questions or submit an issue, if you are unable to run the pipeline locally.

## Running Framework on CTCMS Cluster
Documentation for running pipeline on CTCMS using tmux and the ctcms profile.
CTCMS snakemake cluster deployment profile generated using [Snakemake Slurm profile template](https://github.com/Snakemake-Profiles/slurm).

__Steps__

1. Log into CTCMS using `ssh username@ruth.nist.gov` or `ssh ctcms` if `.ssh/config` is setup to do so (ask Nate O. if you want to do this).
1. Create `tmux` session `tmux new-session -s [session name]`. This will create a detachable terminal session so that the pipeline does not fail if the ssh connection dropped.
1. Switch to appropriate git branch `git checkout [branchname]`, and make sure up to date `git status` and `git pull`. 
1. Activate conda environment for running pipeline `conda activate defrabb`. This environment was create by Nate O. using `mamba env create -f envs/env.yml` from the root directory of this repository.
1. Use complete `config/analyses.tsv` with required run information, update `config/resources.yml` if necessary.
1. Run pipeline using `sh etc/run-analysis-ctcms.sh`


### Notes / Gotchas

* snakemake is executed on headnode as job nodes are not connected to the internet, which is required for conda environments, snakemake wrappers, and downloading resource files.
* Need to first create conda environments using `snakemake --use-conda --conda-create-envs-only`
* Can not define job memory requirement, for jobs with high memory requirements try increasing the number of threads.

## Running Framework on AWS with Tibanna
### Setting up an instance on AWS to run defrabb
Low memory jobs, defined in Snakefile, are run on "local" instance. Higher memory, compute intensive jobs will be run on instances started by Tibanna with appropriate resources. 
1. Set up AWS instance [^2] e.g. i3.large w/ 2vCPU, 15GB and 1TB storage (less storage migth be possible if fewer assemblies will be run)
2. Start instance. The ssh connection command can be found under the "Connect" button in the EC2 console.  Note: you will need to add the path to your `user.pem` file\
`ssh -v -i "~/.ssh/user.pem" ec2-user@##.###.##.##`
3. install dependencies
	- miniconda, find appropriate version of miniconda [here](https://docs.conda.io/en/latest/miniconda.html)\
`wget https://repo.anaconda.com/miniconda/<your_miniconda>.sh`\
`sh <your_miniconda>.sh`\
`echo $HOME` e.g. (/home/ec2-user)\
`export PATH="$HOME/miniconda3/bin:$PATH"`
	- tmux\
`sudo yum install tmux`
	- git\
`sudo yum install git`
	- mamba\
`conda install -c conda-forge mamba`
	- nano (this might already be installed)
4. initialize conda for shell interaction
`conda init bash` (needed to restart shell after this)
5. clone defrabb git repo and switch desired branch\
`git clone https://gitlab.nist.gov/gitlab/njd2/giab-asm-bench-whole-genome.git`\
`git checkout desired-branch`
6. use mamba to set up defrabb environment\
`mamba env create -f envs/env.yml -n defrabb`
7. Add AWS credentials\
get "DEFAULT" credentials ,e.g., `cat ~/.aws/credentials`\
`aws configure` entering the following from your file above when prompted
	- AWS Access Key ID
	- AWS Secret Access Key
	- Default region name 
	- Default output format

### Starting a run on AWS instance that has been configured (see above)
1. Set up `your-directory` using name specific to run, in S3 bucket `giab-tibanna-runs` for run output
2. Setup `config/analyses.tsv`, see `schema/analyses-schema.yml` for description of fields in `analyses.tsv`
3. Setup `config/resources.yml` with specifics for the run(s) outlined in `config/analyses.tsv`.  Make sure to review dipcall and hap.py resources at the bottom of file.
4. Make any necessary changes to `etc/run_on_tibanna_giab.sh`
	- DRYRUN = comment/uncomment, it is suggested you start with a dry run first
	- ANALYSES = add path to file, e.g. `config/myfile.tsv`, if you want to overide use of `anslyses.tsv` defined in `resources.yml` 
	- RUNDIR = `your-directory` in S3 bucket `giab-tibanna-runs` that outputs should go to
	- JOBS and DISKMB = adjust as appropriate for your run requirements
	- set tibanna unicorn  
		If running on AWS from local terminal you can set/use unicorn two ways. Unicorn we are currently using for most work is "tibanna_unicorn_giab_test3"  
		
		**OPTION 1, use if going to consistenly use the same unicorn:**  
		Add the following to your ~ /.bash_profile  
			`export TIBANNA_DEFAULT_STEP_FUNCTION_NAME=<insert unicorn name>`  
		In run_on_tibanna.sh change as follows such that default is used:  
			`PROFILE="default"`  
			`#UNICORN="<insert unicorn name>"`  
			`#export TIBANNA_DEFAULT_STEP_FUNCTION_NAME=${UNICORN}`  
		
		**OPTION 2, set in shell script if you tend to use different unicorns or are running from AWS instance**  
		In run_on_tibanna.sh change as follows such that unicorn is defined:  
			`#PROFILE="default"`  
			`UNICORN="<insert unicorn name>"`  
			`export TIBANNA_DEFAULT_STEP_FUNCTION_NAME=${UNICORN}`  
	
	- review snakemake command and options used by snakemake and tibanna to ensure they are appropriate for your run
5. start tmux session. See [online tmux cheatsheet](https://tmuxcheatsheet.com) for helpful tmux commands [^3]\
`tmux new-session -s my-session-name` multiple people can log in to the session using `tmux a -t my-session-name`
6. activate defrabb environment\
`conda activate defrabb`
7. start run w/ tibanna\
`sh etc/run_on_tibanna_giab.sh`

### Monitoring run on AWS instance
1. Switch from tmux session running script to new tmux session[^3]\
`ctrl + b n` (n=next window) this switches between session windows and can be used to get back to your session.  If timeout occured and you need to log back in to EC2 you can re-attach session after login using `tmux a -t my-session-name` where `my-session-name` is the name of the session you started the run in.
2. Viewing tibanna logs for each job.  This is helpful because job information in the job seesion window goes away as new jobs are started. This will allow you to see if jobs have completed or failed and hopefully gives you information on why they failed. 
	- list the jobs, in this example 10\
`tibanna stat -s tibanna_unicorn_giab_test3 -n 10`\
	- job IDs will be on far left, copy jobID\
`tibanna log -s tibanna_unicorn_giab_test3 -j <paste jobID#>`
3. To view snakemake output while the pipeline is running view snakemake log.  
`cd less .snakemake/log/`\
`less less YYYY-MM-DD<somenumber>.snakemake.log`
4. You can go back to session running script using `ctrl + b n`
5. Monitor memory usage with Tibanna `plot_metrics`.  Note: need newer version of Tibanna than what is available in conda pkg that is installed in defrabb env.  
	a. start local terminal and run the following\
	b. `mamba install pip`\
	c. `pip install -U tibanna`\
	d. `tibanna --version` (currently using tibanna 1.9.1 for this)\
	e. `tibanna plot_metrics -j <insert tibanna jobID#>` this will pull up and html page in browser with metrics

### Notes / Gotchas
- Tibanna and all jobs started with tibanna are run in docker containters.  There is a limit on how many docker containers can be started in a given time. For anonymous login (which we are currently using), limit is 100 pull/6hr, see [Docker pull limit documenation](https://www.docker.com/blog/checking-your-current-docker-pull-rate-limits-and-status/) for mor info on this.  
- If you are downloading numerous large files and run into issues with downloads failing the instance might be out of storage and you might have to increase size of storage through EC2 console after stopping instance. 
- Had a job that failed with no info in tibanna log. The issue was the spot instance was shut down and there wasn't another "reasonable" instance to use.  Users can control how Tibanna uses instances on AWS by adjusting the snakemake command in `etc/run_on_tibanna_giab.sh`. Simply restarting run if it fails migth resolve situation. 
- If hap.py or dipcall fails its worth checking tibanna plot_metrics to make sure they are being given enough memory.  If not adjustments can be made in the `resources.yml`
