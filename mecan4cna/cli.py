import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

import mecan4cna.algorithms as alg
import click 
import mecan4cna.common as comm
from distutils.dir_util import copy_tree
from shutil import copyfile
import subprocess

@click.command()
@click.option('-i', '--input_file', type=click.File('r'), help='The input file.')
@click.option('-o', '--output_path', help='The path for output files.')
@click.option('-p', '--plot', is_flag=True, default=False, help='Whether to show the signal histogram.')
@click.option('-b', '--bins_per_interval', type=click.IntRange(1,None), default=20, help='The number of bins in each copy number interval.')
@click.option('-v', '--intervals', type=click.IntRange(1,None), default=4, help='The number of copy number intervals.')
@click.option('--demo', is_flag=True, help='Copy example files and run a demo script in the current directry.')
@click.option('-pt', '--peak_thresh', type=click.IntRange(0,None), default=1000, help='The minimum probes of a peak.')
@click.option('-st', '--segment_thresh', type=click.IntRange(0,None), default=3, help='The minimum probes of a segment.')
@click.option('--model_steps', type=click.IntRange(0,None), default=2, help='The incremental step size in modeling.')
@click.option('--mpd_coef', type=click.FLOAT, default=0.1, help='Minimun Peak Distance coefficient in peak detection.')
@click.option('--max_level_distance', type=click.FLOAT, default=1.3, help='The maximum value of level distance.')
@click.option('--min_level_distance', type=click.FLOAT, default=0.3, help='The minimum value of level distance.')
@click.option('--min_model_score', type=click.IntRange(0,None), default=9, help='The minimum value of model score.')
@click.option('--info_lost_ratio_thresh', type=click.FLOAT, default=0.3, help='The threshold of information lost ratio.')
@click.option('--info_lost_range_low', type=click.FLOAT, default=0.2, help='The low end of information lost range.')
@click.option('--info_lost_range_high', type=click.FLOAT, default=0.8, help='The high end of information lost range.')
def cli(input_file,output_path,plot,bins_per_interval,intervals,peak_thresh,segment_thresh,model_steps,mpd_coef,max_level_distance,min_level_distance,
        min_model_score, info_lost_ratio_thresh, info_lost_range_low,info_lost_range_high, demo ):
    

    if demo:
        example_path = 'examples'
        script_name = 'run_mecan_example.sh'
        examples_path = os.path.join(os.path.dirname(__file__), example_path)
        script_path = os.path.join(os.path.dirname(__file__), script_name)
        copy_tree(examples_path, example_path, preserve_mode=0)
        copyfile(script_path , script_name)
        if os.path.isdir(example_path):
            print('Copied examples to {}'.format(os.path.abspath(example_path)))
            print('Copied script to {}'.format(os.path.abspath('.')))

        subprocess.run(['sh',script_name], stdout=sys.stdout, stderr=sys.stderr)
        exit()
    

    if input_file:
        segments = comm.file2list(input_file)
    else:
        exit("Specify an input file.")

    if output_path:
        os.makedirs(output_path, exist_ok=True)
    else:
        exit("Specify an output path.")



    r = alg.mecan(bins_per_interval=bins_per_interval, intervals=intervals, interval_step=model_steps, minprobes=peak_thresh,
        probe_thresh=segment_thresh, mpd_coef=mpd_coef, outpath=output_path, plot=plot, thresh_min=min_level_distance, thresh_max=max_level_distance, 
        minlevel=min_model_score, neglect_level_cap_high= info_lost_range_high, neglect_level_cap_low=info_lost_range_low, neglect_thresh=info_lost_ratio_thresh)

    res = r.run(segments)

def main():
    try:
        cli()
    except Exception as e:
        print(e)

if __name__ == '__main__':
    main()