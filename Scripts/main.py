import subprocess
import os

def submit_job(script_path):
    """Submit bash scripts as jobs"""
    try:
        subprocess.run(['sbatch', script_path], check=True)
        print(f"Successfully submitted {script_path} using sbatch")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while submitting {script_path}: {e}")

def main():
    # Paths to the scripts
    script1 = os.path.join('CountMatrix', 'run_count.sh')
    script2 = os.path.join('CellBender', 'run_cellbender.sh')

    # Submit scripts in order
    submit_job(script1)
    submit_job(script2)

if __name__ == "__main__":
    main()
