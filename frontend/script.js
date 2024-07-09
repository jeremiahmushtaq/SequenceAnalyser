function uploadFiles(inputId, formId) {
    const fileInput = document.getElementById(inputId);
    const tickMark = document.getElementById(`${inputId.replace('-upload', '')}-tick`);
    const form = document.getElementById(formId);
    
    fileInput.click();
    
    fileInput.onchange = () => {
        if (fileInput.files.length > 0) {
            tickMark.style.display = 'inline';

            const formData = new FormData();
            for (let i = 0; i < fileInput.files.length; i++) {
                formData.append(fileInput.name, fileInput.files[i]);
            }

            fetch(`/upload_${inputId.replace('-upload', '')}`, {
                method: 'POST',
                body: formData
            }).then(response => response.json())
              .then(data => console.log(data))
              .catch(error => console.error('Error:', error));
        } else {
            tickMark.style.display = 'none';
        }
    };
}