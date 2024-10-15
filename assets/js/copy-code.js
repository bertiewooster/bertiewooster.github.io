document.addEventListener("DOMContentLoaded", function() {
    const codeBlocks = document.querySelectorAll("pre code");
  
    codeBlocks.forEach(codeBlock => {
      const button = document.createElement("button");
      button.className = "copy-btn";
      button.textContent = "Copy";
  
      const pre = codeBlock.parentElement;
      pre.parentElement.insertBefore(button, pre);
  
      button.addEventListener("click", function() {
        const textarea = document.createElement("textarea");
        textarea.value = codeBlock.textContent;
        document.body.appendChild(textarea);
        textarea.select();
        document.execCommand('copy');
        document.body.removeChild(textarea);
  
        button.textContent = "Copied!";
        setTimeout(() => {
          button.textContent = "Copy";
        }, 2000);
      });
    });
  });
  