The step-by-step workflow to convert your list of compound names into CIDs or structure data using only the web interface is detailed below:

Step 1: Format your text file
Make sure your text file containing the compound names is formatted properly:

Open your list in a basic text editor (like Notepad or TextEdit).

Ensure there is only one compound name per line.

Do not include any headers (e.g., do not write "Compound Name" on the first line).

Save it as a plain text file (.txt).
(Note: If you have a column of names in Excel, you can just copy that column and paste it directly into the web form in the next step).

Step 2: Use the PubChem Identifier Exchange Service
Go to the PubChem Identifier Exchange Service page.

Under 1) Select an input format and give the ID list:

Change the dropdown from Registry ID to Synonyms (this is what PubChem calls chemical/compound names).

Either click Browse... to upload your .txt file, or simply copy and paste your list of names into the large text box.

Under 2) Select an operation type:

Leave this as Same CID. This tells the system to find the PubChem Compound ID that matches your name.

Under 3) Select an output type:

Select CIDs if you just want a list of IDs to feed into other tools.

Tip: If your end goal is getting the structure files (like SMILES or InChIKeys) for these names, you can actually select SMILES or InChIKeys right here to skip an extra step!

Under 4) Select an output method:

Select Two column file showing input-output correspondence. This is highly recommended because it will give you a file showing [Your Input Name] -> [PubChem CID]. If a name wasn't found, it will let you know, so your list doesn't get mixed up.

Scroll down to the bottom and click Submit Job.

Step 3: Download the Results
After clicking submit, the page will refresh to a waiting screen. Depending on how long your list is (it handles up to 500,000 names), it will take anywhere from a few seconds to a few minutes.

Once finished, a download link to a text file will appear.

Step 4: (Optional) Feed the CIDs into the Structure Download Service
If you selected "CIDs" in Step 3 but your ultimate goal is to get 2D/3D structure files (like .sdf files) for your molecular modeling software:

Open the downloaded text file and copy the column containing the new CIDs.

Go back to the PubChem Structure Download Service.

Paste the CIDs there, choose SDF (or your format of choice), and hit download.