import os
import sys
import time
import logging
from pathlib import Path

# Registra el tiempo de inicio
start_time = time.time()

import frustraevo.Functions as Functions
import argparse

logger = logging.getLogger(__name__)

def main(argv):

    parser = argparse.ArgumentParser(description='Calculation of the frustration logo.')
    parser.add_argument("--JobId", help="The name of the job")
    parser.add_argument("-o", help="path to the output directory. (should not exist before running the script)")
    #parser.add_argument("--RPath", default='N', help="Path to R script files (Default: Scripts)")
    parser.add_argument("--fasta", help="Name of the fasta")
    parser.add_argument("--ref", default='None', help="Id of the reference protein of your logo (Default: None)")
    parser.add_argument("--pdb_db", default='None', help="Path to the PDBs folder (Default: None)")
    parser.add_argument("--cmaps", default=False, action="store_true", help="Set to activate contact maps calculation")
    

    #How to run the pipeline in linux terminal:
    #python run_logo.py JobId RPath fasta list_file ref pdb_db 

    args = parser.parse_args(argv)
    list_file=''

    out_dir = Path(args.o)

    RPath = Path(__file__).parent / "Scripts"

    if out_dir.exists():
        raise ValueError(f"Output directory {out_dir} already exists. Please remove it before running the script.")
    
    out_dir.mkdir()

    logger.basicFonfig(filename=out_dir/"log.txt")
        #the parameter Scripts is the path to the .r and .py files for plots
    logger.info('Copying files for Frustration Logo')
    Functions.setup_outdir(out_dir)
    
    list_file =  out_dir / "prot_names.txt"
    Functions.pdb_list(args.fasta, list_file)
    
    logger.info('Changes in MSA Frustration Logo')
    Functions.changes(out_dir, args.fasta)
    
    logger.info('Running Checks in sequence')
    Functions.checks_seq(list_file, out_dir, args.pdb_db)
    
    path_file='FrustraEvo_'+args.JobId+'/ErrorSeq.log'
    logger.info(' Reading and Preparing Files...')
    
    if not os.path.exists(path_file):
        logger.info(' Calculating Single Residue Frustration Index...')
        Functions.FrustraPDB(list_file,args.JobId,args.pdb_db)
        logger.info(' Running Internal Checks...')

        Functions.checks(args.JobId)
        logger.info('Preparing MSA Files...')
        
        Functions.prepare_file(args.JobId,args.ref)
        Functions.FinalAlign(args.JobId)
        logger.info('Running Equivalences...')
        
        Functions.Equivalences(args.JobId)
        logger.info("Calculating SeqIC and FrustIC...")
        Functions.FastaMod(args.JobId)
        logger.info("Running Checks")
        Functions.LogoCheck(args.JobId)
        logger.info('Generating Sequence and Frustration Logos Plots...')
        Functions.plots_logo(args.JobId,args.ref,args.RPath) 
        
        if args.cmaps:
            logger.info('Calculating Mutational Frustration Index and Contact Maps...')
            Functions.CMaps_Mutational(args.JobId,args.RPath,args.ref)#Genera los mapas de contacto para el indice mutational
            logger.info("Calculating Configurational Frustration Index and Contact Maps...")
            Functions.CMaps_Configurational(args.JobId,args.RPath,args.ref)#Genera los mapas de contacto para el indice configurational
        Functions.clean_files(args.JobId,args.RPath,args.ref)    
        end_time = time.time()
        elapsed_time = end_time - start_time
        logger.info(f"Elapsed time: {elapsed_time}")
        logger.info(f"Generating Output Files...")
        Functions.VScript(args.JobId,elapsed_time)
        os.system('cd '+args.RPath+';python3 setup_render.py '+args.JobId+' SingleRes')
        os.system('cd '+args.RPath+';python3 Seq_IC.py '+args.JobId)
    #contact_maps.py	
        os.system('cd '+args.RPath+';python3 contact_maps.py '+args.JobId+' '+args.ref+' IC_SingleRes_'+args.JobId+' IC_Mut_'+args.JobId+' IC_Conf_'+args.JobId)

        logger.info("Job finished")

if __name__ == "__main__":
        main(sys.argv[1:])

def _entrypoint():
        main(sys.argv[1:])