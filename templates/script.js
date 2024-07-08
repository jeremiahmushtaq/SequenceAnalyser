// script.js

function uploadFile(inputId) {
    const fileInput = document.getElementById(inputId);
    const tickMark = document.getElementById(`${inputId.replace('-upload', '')}-tick`);
    
    // Simulate a click on the hidden file input
    fileInput.click();
    
    fileInput.onchange = () => {
        if (fileInput.files.length > 0) {
            // Show the green tick mark when a file is selected
            tickMark.style.display = 'inline';
        } else {
            // Hide the tick mark if no file is selected
            tickMark.style.display = 'none';
        }
    };
}