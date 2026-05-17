import requests
import gzip
import shutil
import os
import concurrent.futures
import time

def download_protein(protein_code, session):
    url = f"https://files.rcsb.org/download/{protein_code}.pdb1.gz"
    cwd = os.getcwd()

    # Tạo thư mục "protein" nếu nó chưa tồn tại
    protein_dir = os.path.join(cwd, "protein")
    os.makedirs(protein_dir, exist_ok=True)

    try:
        # Add sleep time to prevent rate limiting or overwhelming the server
        time.sleep(1)

        # Send a GET request to the URL using the provided session
        response = session.get(url)

        # Check if the request was successful
        if response.status_code == 200:
            # Define the filenames with path (now includes the "protein" directory)
            gz_filename = os.path.join(protein_dir, f"{protein_code}.pdb1.gz")
            pdb_filename = os.path.join(protein_dir, f"{protein_code}.pdb")

            # Save the content to a .gz file
            with open(gz_filename, 'wb') as file:
                file.write(response.content)

            print(f"Protein {protein_code} downloaded successfully as {gz_filename}.")

            # Unzip the .gz file
            with gzip.open(gz_filename, 'rb') as f_in:
                with open(pdb_filename, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            print(f"Protein {protein_code} unzipped successfully as {pdb_filename}.")

            # Optionally, delete the .gz file
            os.remove(gz_filename)

        else:
            print(f"Failed to download protein {protein_code}. HTTP Status Code: {response.status_code}")

    except Exception as e:
        print(f"An error occurred: {e}")

def main():
    cwd = os.getcwd()
    protein_file_path = os.path.join(cwd, "protein.txt")
    print(f"Current working directory: {cwd}")
    print(f"Protein file path: {protein_file_path}")

    try:
        with open(protein_file_path, 'r') as file:
            protein_codes = [line.strip() for line in file]
            print(f"Protein codes: {protein_codes}")

        # Create a requests session
        with requests.Session() as session:
            # Use ThreadPoolExecutor for multithreading
            with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
                # Submit download tasks to the executor
                futures = [executor.submit(download_protein, protein_code, session) for protein_code in protein_codes]
                # Wait for all tasks to complete and handle potential exceptions
                for future in concurrent.futures.as_completed(futures):
                    try:
                        future.result()
                    except Exception as e:
                        print(f"Error during download: {e}")

    except FileNotFoundError:
        print(f"Error: File 'protein.txt' not found at {protein_file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()