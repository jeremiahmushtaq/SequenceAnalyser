from flask import Flask, request, render_template, redirect
from werkzeug.utils import secure_filename
import os
from functions.primers_and_melting_temperature import sequencePCRtemp, BadSequenceException
from functions.translation_and_reading_frames import maximalORFa, OpenReadingFrameException
from functions.pairwise_distance_matrix import DistanceMatrix, DimensionalityException
from werkzeug.exceptions import RequestEntityTooLarge

app = Flask(__name__)
app.config['UPLOAD_DIRECTORY'] = 'uploads'  # Directory to store uploaded files
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # Maximum file size limit set to 16 MB
app.config['ALLOWED_EXTENSIONS'] = ['.txt', '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn']  # Allowed file extensions

@app.route('/')
def index():
    """Render the main upload page."""
    return render_template('index.html')

@app.route('/upload/primer', methods=['POST'])
def Calculate_Prime_Melting_Temperature():
    """
    Handle the file upload, ensure the file is valid, process it,
    render the results, and delete the file afterward.
    """
    try:
        # Retrieve the uploaded file from the form
        file = request.files['file']
        
        # Extract the file extension to validate the file type
        extension = os.path.splitext(file.filename)[1]

        # Check if the file exists and has an allowed extension
        if file:
            if extension not in app.config['ALLOWED_EXTENSIONS']:
                # Return error for invalid file types
                return 'File is not a TXT or FASTA file.'  

            # Securely save the file to the designated uploads directory
            file_path = os.path.join(
                app.config['UPLOAD_DIRECTORY'],
                secure_filename(file.filename)
            )
            file.save(file_path)

            try:
                # Process the uploaded file using the custom function
                result = sequencePCRtemp(file_path)

            except BadSequenceException as e:
                # Handle the specific exception raised by your function
                return f"{str(e)}"

            # Render the results in the results.html template
            return render_template('results_primer.html', result=result)

    except RequestEntityTooLarge:
        # Handle file size errors
        return 'File is larger than 16MB limit.'
    
    finally:
        # Ensure that the file is deleted, whether or not an exception occurred
        if file_path and os.path.exists(file_path):
            try:
                os.remove(file_path)
                print(f"Deleted file: {file_path}")  # Log the file deletion
            except OSError as e:
                print(f"Error deleting file {file_path}: {e}")  # Log any errors in file deletion

    # Redirect to the home page if no file is uploaded
    return redirect("/")

@app.route('/upload/orf', methods=['POST'])
def Longest_Open_Reading_Frame():
    """
    Handle the file upload, ensure the file is valid, process it,
    render the results, and delete the file afterward.
    """
    try:
        # Retrieve the uploaded file from the form
        file = request.files['file']
        
        # Extract the file extension to validate the file type
        extension = os.path.splitext(file.filename)[1]

        # Check if the file exists and has an allowed extension
        if file:
            if extension not in app.config['ALLOWED_EXTENSIONS']:
                # Return error for invalid file types
                return 'File is not a TXT or FASTA file.'  

            # Securely save the file to the designated uploads directory
            file_path = os.path.join(
                app.config['UPLOAD_DIRECTORY'],
                secure_filename(file.filename)
            )
            file.save(file_path)

            try:
                # Process the uploaded file using the custom function
                result = maximalORFa(file_path)

            except BadSequenceException as e:
                # Handle the specific exception raised by your function
                return f"{str(e)}"
            
            except OpenReadingFrameException as e:
                # Handle the specific exception raised by your function
                return f"{str(e)}"

            # Render the results in the results.html template
            return render_template('results_orf.html', result=result)

    except RequestEntityTooLarge:
        # Handle file size errors
        return 'File is larger than 16MB limit.'
    
    finally:
        # Ensure that the file is deleted, whether or not an exception occurred
        if file_path and os.path.exists(file_path):
            try:
                os.remove(file_path)
                print(f"Deleted file: {file_path}")  # Log the file deletion
            except OSError as e:
                print(f"Error deleting file {file_path}: {e}")  # Log any errors in file deletion

    # Redirect to the home page if no file is uploaded
    return redirect("/")

@app.route('/upload/pdm', methods=['POST'])
def Pairwise_Distance_Matrix():
    """
    Handle the file upload, ensure the file is valid, process it,
    render the results, and delete the files afterward.
    """
    file_paths = []  # List to store paths of saved files

    try:
        # Retrieve the uploaded files from the form
        files = request.files.getlist('file')

        if not files:
            return 'No files uploaded.'

        for file in files:
            # Extract the file extension to validate the file type
            extension = os.path.splitext(file.filename)[1]

            # Check if the file exists and has an allowed extension
            if file:
                if extension not in app.config['ALLOWED_EXTENSIONS']:
                    # Return error for invalid file types
                    return 'One of the files is not a valid TXT or FASTA file.'  

                # Securely save the file to the designated uploads directory
                file_path = os.path.join(
                    app.config['UPLOAD_DIRECTORY'],
                    secure_filename(file.filename)
                )
                file.save(file_path)
                file_paths.append(file_path)  # Append saved file path to the list

        # List the saved files (full paths) to pass to the DistanceMatrix function
        user_files = file_paths

        try:
            # Process the uploaded files using the custom DistanceMatrix function
            result = DistanceMatrix(user_files)

        except DimensionalityException as e:
            # Handle the specific exception raised by your function
            return f"Error: {str(e)}"

        # Render the results in the results_pdm.html template
        return render_template('results_pdm.html', result=result)

    except RequestEntityTooLarge:
        # Handle file size errors
        return 'One or more files are larger than the 16MB limit.'

    finally:
        # Ensure that the uploaded files are deleted after processing
        for file_path in file_paths:
            if os.path.exists(file_path):
                try:
                    os.remove(file_path)
                    print(f"Deleted file: {file_path}")
                except OSError as e:
                    print(f"Error deleting file {file_path}: {e}")

    # Redirect to the home page if no file is uploaded
    return redirect("/")

if __name__ == '__main__':
    app.run(debug=True)