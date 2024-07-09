from flask import Flask, request, jsonify, send_from_directory
import os

app = Flask(__name__, static_folder='frontend')
UPLOAD_FOLDER = 'uploads'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route('/')
def serve_index():
    return send_from_directory(app.static_folder, 'index.html')

@app.route('/<path:path>')
def serve_static(path):
    return send_from_directory(app.static_folder, path)

@app.route('/upload_primer', methods=['POST'])
def upload_primer():
    return handle_file_upload('primer-upload')

@app.route('/upload_amino', methods=['POST'])
def upload_amino():
    return handle_file_upload('amino-upload')

@app.route('/upload_orf', methods=['POST'])
def upload_orf():
    return handle_file_upload('orf-upload')

def handle_file_upload(file_key):
    if file_key not in request.files:
        return "No file part", 400
    files = request.files.getlist(file_key)
    if not files:
        return "No selected files", 400

    uploaded_files = []
    for file in files:
        if file and file.filename != '':
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
            file.save(filepath)
            uploaded_files.append(file.filename)

    return jsonify({"message": "Files uploaded successfully", "filenames": uploaded_files})

if __name__ == '__main__':
    app.run(debug=True)