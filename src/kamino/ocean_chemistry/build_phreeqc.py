import subprocess
import shutil
import tarfile
import sys
from pathlib import Path

def build():
    # 1. Define Directories
    # BASE_DIR is src/planet_model/ocean/
    BASE_DIR = Path(__file__).parent.absolute()
    
    # We will do the work in a temporary folder
    BUILD_DIR = BASE_DIR / "build_temp"
    
    # We will install the final result here: src/planet_model/ocean/phreeqc_bin
    INSTALL_DIR = BASE_DIR / "phreeqc_bin"
    
    # The source tarball (Must be in the same folder as this script)
    TARBALL_NAME = "phreeqc-3.8.6-17100.tar.gz"
    TARBALL_PATH = BASE_DIR / "phreeqc" / TARBALL_NAME

    # 2. Safety Checks
    if not TARBALL_PATH.exists():
        print(f"Error: Could not find {TARBALL_NAME} in {BASE_DIR}")
        print("Please download the source tarball and place it in this folder.")
        sys.exit(1)

    # Clean up previous builds to ensure a fresh start
    if BUILD_DIR.exists():
        shutil.rmtree(BUILD_DIR)
    if INSTALL_DIR.exists():
        shutil.rmtree(INSTALL_DIR)

    BUILD_DIR.mkdir()
    
    try:
        print("--- Step 1: Extracting Source ---")
        with tarfile.open(TARBALL_PATH, "r:gz") as tar:
            tar.extractall(path=BUILD_DIR)
        
        # Find the extracted folder name (e.g., phreeqc-3.8.6-17100)
        # We look for the first directory inside build_temp
        extracted_folder_name = [f.name for f in BUILD_DIR.iterdir() if f.is_dir()][0]
        SOURCE_DIR = BUILD_DIR / extracted_folder_name
        
        # Create the 'Release' folder exactly like the bash script did
        RELEASE_DIR = SOURCE_DIR / "Release"
        RELEASE_DIR.mkdir()

        print("--- Step 2: Configuring ---")
        # equivalent to: ../configure --prefix=$INSTALLDIR
        # We use absolute paths to avoid confusion
        configure_script = SOURCE_DIR / "configure"
        
        subprocess.run(
            [str(configure_script), f"--prefix={INSTALL_DIR}"],
            cwd=RELEASE_DIR,
            check=True
        )

        print("--- Step 3: Compiling (Make) ---")
        # Run make inside the Release folder
        # '-j' uses multiple cores for faster compilation
        subprocess.run(["make", "-j"], cwd=RELEASE_DIR, check=True)

        print("--- Step 4: Installing ---")
        subprocess.run(["make", "install"], cwd=RELEASE_DIR, check=True)

        print("\nSUCCESS!")
        print(f"PHREEQC installed to: {INSTALL_DIR}/bin/phreeqc")

    except subprocess.CalledProcessError as e:
        print(f"\nError during compilation step: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\nUnexpected error: {e}")
        sys.exit(1)
    finally:
        # Optional: Cleanup the messy build folder, keep the install folder
        if BUILD_DIR.exists():
            print("Cleaning up temporary build files...")
            shutil.rmtree(BUILD_DIR)

if __name__ == "__main__":
    build()