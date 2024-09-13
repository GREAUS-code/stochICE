import os
import shutil
import subprocess

def copy_folders_and_modify_script(src_folder, script_path, n):
    # Get the parent directory and the folder name from the specified path
    parent_dir = os.path.dirname(src_folder)
    folder_name = os.path.basename(src_folder)
    script_dir = os.path.dirname(script_path)
    script_name = os.path.basename(script_path)
    scripts = []

    for i in range(1, n + 1):
        # Define the new folder name with '_n' appended
        new_folder_name = f"{folder_name}_{i}"
        new_folder_path = os.path.join(parent_dir, new_folder_name)
        new_script_name = f"{os.path.splitext(script_name)[0]}_{i}.py"

        # Copy the folder to the same level with the new name
        shutil.copytree(src_folder, new_folder_path)
        print(f"Copied folder '{src_folder}' to '{new_folder_path}'")

        # Read the script and modify the `path` variable
        with open(script_path, 'r') as script_file:
            script_content = script_file.read()

        # Modify the path to reflect the new folder name
        modified_script_content = script_content.replace(
            "r'C:\\Users\\Jason\\Desktop\\stochICE\\examples\\Chateauguay'", 
            f"r'{new_folder_path}'"
        )

        # Write the modified script to the same directory as the original script
        new_script_path = os.path.join(script_dir, new_script_name)
        with open(new_script_path, 'w') as new_script_file:
            new_script_file.write(modified_script_content)
        print(f"Copied and modified script '{script_name}' to '{new_script_path}'")

        # Add the new script name to the list
        scripts.append(new_script_path)

    return scripts

def run_scripts_in_parallel(scripts):
    processes = []
    
    # Launch all the modified scripts in parallel
    for script in scripts:
        print(f"Launching {script}...")
        process = subprocess.Popen(["python", script])
        processes.append(process)
    
    # Wait for all processes to finish
    for process in processes:
        process.wait()

def copy_generated_data_back(src_folder, n):
    # Define subfolders in MonteCarlo directory
    monte_carlo_subfolders = ['FlowFiles', 'GeoFiles', 'SimulationTifs']
    
    # Loop through each generated folder
    for i in range(1, n + 1):
        # Define the new folder name with '_n' appended
        new_folder_name = f"{os.path.basename(src_folder)}_{i}"
        new_folder_path = os.path.join(os.path.dirname(src_folder), new_folder_name)
        
        for subfolder in monte_carlo_subfolders:
            src_subfolder_path = os.path.join(new_folder_path, 'MonteCarlo', subfolder)
            dest_subfolder_path = os.path.join(src_folder, 'MonteCarlo', subfolder)
            
            if os.path.exists(src_subfolder_path):
                # Copy the contents of the generated subfolder back to the original subfolder
                for item in os.listdir(src_subfolder_path):
                    src_item = os.path.join(src_subfolder_path, item)
                    dest_item = os.path.join(dest_subfolder_path, item)
                    
                    if os.path.isfile(src_item):
                        shutil.copy2(src_item, dest_item)
                    elif os.path.isdir(src_item):
                        # If the item is a directory, use copytree
                        if os.path.exists(dest_item):
                            shutil.rmtree(dest_item)  # Remove existing folder first
                        shutil.copytree(src_item, dest_item)
                
                print(f"Copied contents from '{src_subfolder_path}' to '{dest_subfolder_path}'")
            else:
                print(f"Subfolder '{src_subfolder_path}' does not exist")

if __name__ == "__main__":
    # Specify the source folder, the script path, and number of copies
    src_folder = r"C:\Users\Jason\Desktop\stochICE\examples\Chateauguay"  # Absolute path of the source folder
    script_path = r"C:\Users\Jason\Desktop\stochICE\examples\stoch_Chateauguay.py"  # Path to the original Python script
    n = 5  # Number of copies

    # Copy folders and modify scripts
    print("Copying folders and modifying scripts...")
    scripts_list = copy_folders_and_modify_script(src_folder, script_path, n)

    # Ensure everything is copied and modified before launching the scripts
    if scripts_list:
        print("\nAll scripts copied and modified successfully. Starting parallel execution...\n")
        run_scripts_in_parallel(scripts_list)
        
        print("\nParallel execution finished. Copying generated data back to original folder...\n")
        copy_generated_data_back(src_folder, n)
        
        print("\nData copied back to original folder successfully.")
    else:
        print("No scripts were copied or modified. Aborting execution.")
