document.addEventListener('DOMContentLoaded', function() {
    const generateBtn = document.getElementById('generate-btn');
    const artImage = document.getElementById('art-image');
    const seedValue = document.getElementById('seed-value');
    const loading = document.getElementById('loading');
    const galleryItems = document.getElementById('gallery-items');
    const downloadBtn = document.getElementById('download-btn');
    const animateBtn = document.getElementById('animate-btn');
    
    // Generate initial artwork
    generateArt();
    
    // Event listener for generate button
    generateBtn.addEventListener('click', generateArt);
    
    // Function to generate art
    function generateArt() {
        // Show loading spinner
        loading.style.display = 'block';
        artImage.style.opacity = '0.3';
        
        // Get selected options
        const style = document.getElementById('style').value;
        const colormap = document.getElementById('colormap').value;
        const seed = document.getElementById('seed').value;
        
        // Build URL with parameters
        let url = `/generate?style=${style}&cmap=${colormap}`;
        if (seed) {
            url += `&seed=${seed}`;
        }
        
        // Fetch the generated art
        fetch(url)
            .then(response => {
                if (!response.ok) {
                    throw new Error(`HTTP error! Status: ${response.status}`);
                }
                return response.json();
            })
            .then(data => {
                // Update the image
                artImage.src = `data:image/png;base64,${data.image}`;
                seedValue.textContent = data.seed;
                
                // Add to gallery
                addToGallery(data.image, style, colormap, data.seed);
                
                // Hide loading spinner
                loading.style.display = 'none';
                artImage.style.opacity = '1';
            })
            .catch(error => {
                console.error('Error generating art:', error);
                loading.style.display = 'none';
                artImage.style.opacity = '1';
                alert('Error generating artwork. Please try again.');
            });
    }
    
    // Function to add art to gallery
    function addToGallery(imageData, style, colormap, seed) {
        // Create gallery item
        const galleryItem = document.createElement('div');
        galleryItem.className = 'gallery-item';
        
        // Map style names to more readable names
        const styleNames = {
            'spiral': 'Spiral Waves',
            'fractal': 'Fractal Dream',
            'cosmic': 'Cosmic Flow',
            'harmonic': 'Harmonic Structure'
        };
        
        // Map colormap names to more readable names
        const colormapNames = {
            'cosmic': 'Cosmic',
            'sunset': 'Sunset',
            'ocean': 'Ocean',
            'forest': 'Forest'
        };
        
        // Create gallery item HTML
        galleryItem.innerHTML = `
            <img src="data:image/png;base64,${imageData}" alt="Generated Art">
            <div class="gallery-item-info">
                <p><strong>Style:</strong> ${styleNames[style]}</p>
                <p><strong>Colors:</strong> ${colormapNames[colormap]}</p>
                <p><strong>Seed:</strong> ${seed}</p>
                <button class="btn regenerate-btn" data-seed="${seed}" data-style="${style}" data-colormap="${colormap}">Regenerate</button>
            </div>
        `;
        
        // Add to gallery
        galleryItems.prepend(galleryItem);
        
        // Add event listener to regenerate button
        const regenerateBtn = galleryItem.querySelector('.regenerate-btn');
        regenerateBtn.addEventListener('click', function() {
            const seed = this.getAttribute('data-seed');
            const style = this.getAttribute('data-style');
            const colormap = this.getAttribute('data-colormap');
            
            // Set form values
            document.getElementById('style').value = style;
            document.getElementById('colormap').value = colormap;
            document.getElementById('seed').value = seed;
            
            // Generate art with these settings
            generateArt();
            
            // Scroll to the art display
            document.querySelector('.art-container').scrollIntoView({ behavior: 'smooth' });
        });
    }
    
    // Add animation effects for smoother user experience
    artImage.addEventListener('load', function() {
        this.classList.add('fade-in');
        setTimeout(() => {
            this.classList.remove('fade-in');
        }, 500);
    });
    
    // Event listener for download button
    downloadBtn.addEventListener('click', function() {
        // Get the current image
        const imageUrl = artImage.src;
        
        // Send to server for download
        fetch('/download', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({ image: imageUrl }),
        })
        .then(response => {
            if (!response.ok) {
                throw new Error(`HTTP error! Status: ${response.status}`);
            }
            return response.blob();
        })
        .then(blob => {
            // Create a download link
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.style.display = 'none';
            a.href = url;
            a.download = `generative_art_${seedValue.textContent}.png`;
            document.body.appendChild(a);
            a.click();
            window.URL.revokeObjectURL(url);
        })
        .catch(error => {
            console.error('Error downloading image:', error);
            alert('Error downloading image. Please try again.');
        });
    });
    
    // Event listener for animate button
    animateBtn.addEventListener('click', function() {
        // Show loading spinner
        loading.style.display = 'block';
        
        // Get selected options
        const style = document.getElementById('style').value;
        const colormap = document.getElementById('colormap').value;
        
        // Fetch animation frames
        fetch(`/animate?style=${style}&cmap=${colormap}&frames=30`)
            .then(response => {
                if (!response.ok) {
                    throw new Error(`HTTP error! Status: ${response.status}`);
                }
                return response.json();
            })
            .then(data => {
                // Hide loading spinner
                loading.style.display = 'none';
                
                // Create animation container
                const animationContainer = document.createElement('div');
                animationContainer.className = 'animation-container';
                
                // Create animation canvas (img element)
                const canvas = document.createElement('img');
                canvas.className = 'animation-canvas';
                canvas.src = `data:image/png;base64,${data.frames[0]}`;
                
                // Create controls
                const controls = document.createElement('div');
                controls.className = 'animation-controls';
                
                // Close button
                const closeBtn = document.createElement('button');
                closeBtn.className = 'btn';
                closeBtn.innerHTML = '<i class="fas fa-times"></i> Close';
                closeBtn.addEventListener('click', function() {
                    document.body.removeChild(animationContainer);
                });
                
                // Save animation button
                const saveBtn = document.createElement('button');
                saveBtn.className = 'btn';
                saveBtn.innerHTML = '<i class="fas fa-save"></i> Save Frame';
                saveBtn.addEventListener('click', function() {
                    // Get current frame
                    const currentFrame = canvas.src;
                    
                    // Send to server for download
                    fetch('/download', {
                        method: 'POST',
                        headers: {
                            'Content-Type': 'application/json',
                        },
                        body: JSON.stringify({ image: currentFrame }),
                    })
                    .then(response => response.blob())
                    .then(blob => {
                        // Create a download link
                        const url = window.URL.createObjectURL(blob);
                        const a = document.createElement('a');
                        a.style.display = 'none';
                        a.href = url;
                        a.download = `generative_art_animation_${data.base_seed}.png`;
                        document.body.appendChild(a);
                        a.click();
                        window.URL.revokeObjectURL(url);
                    })
                    .catch(error => console.error('Error downloading image:', error));
                });
                
                // Play/pause button
                const playPauseBtn = document.createElement('button');
                playPauseBtn.className = 'btn';
                playPauseBtn.innerHTML = '<i class="fas fa-pause"></i> Pause';
                let isPlaying = true;
                playPauseBtn.addEventListener('click', function() {
                    if (isPlaying) {
                        clearInterval(animationInterval);
                        playPauseBtn.innerHTML = '<i class="fas fa-play"></i> Play';
                    } else {
                        startAnimation();
                        playPauseBtn.innerHTML = '<i class="fas fa-pause"></i> Pause';
                    }
                    isPlaying = !isPlaying;
                });
                
                // Add elements to container
                controls.appendChild(playPauseBtn);
                controls.appendChild(saveBtn);
                controls.appendChild(closeBtn);
                animationContainer.appendChild(canvas);
                animationContainer.appendChild(controls);
                document.body.appendChild(animationContainer);
                
                // Start animation
                let frameIndex = 0;
                let animationInterval;
                
                function startAnimation() {
                    animationInterval = setInterval(() => {
                        frameIndex = (frameIndex + 1) % data.frames.length;
                        canvas.src = `data:image/png;base64,${data.frames[frameIndex]}`;
                    }, 100);
                }
                
                startAnimation();
                
                // Stop animation when container is closed
                closeBtn.addEventListener('click', function() {
                    clearInterval(animationInterval);
                });
            })
            .catch(error => {
                console.error('Error generating animation:', error);
                loading.style.display = 'none';
                alert('Error generating animation. Please try again.');
            });
    });
    
    // Add keyboard shortcuts
    document.addEventListener('keydown', function(event) {
        // Generate new art with spacebar
        if (event.code === 'Space' && !event.target.matches('input, select, textarea')) {
            event.preventDefault();
            generateBtn.click();
        }
        
        // Download with Ctrl+S
        if (event.code === 'KeyS' && (event.ctrlKey || event.metaKey) && !event.target.matches('input, select, textarea')) {
            event.preventDefault();
            downloadBtn.click();
        }
        
        // Animate with Ctrl+A
        if (event.code === 'KeyA' && (event.ctrlKey || event.metaKey) && !event.target.matches('input, select, textarea')) {
            event.preventDefault();
            animateBtn.click();
        }
    });
    
    // Add theme toggle functionality
    const themeToggle = document.createElement('button');
    themeToggle.className = 'theme-toggle';
    themeToggle.innerHTML = '<i class="fas fa-moon"></i>';
    themeToggle.title = 'Toggle Light/Dark Mode';
    document.querySelector('.container').appendChild(themeToggle);
    
    let isDarkMode = true;
    themeToggle.addEventListener('click', function() {
        document.body.classList.toggle('light-mode');
        isDarkMode = !isDarkMode;
        themeToggle.innerHTML = isDarkMode ? '<i class="fas fa-moon"></i>' : '<i class="fas fa-sun"></i>';
    });
    
    // Add responsive gallery view toggle
    const galleryViewToggle = document.createElement('button');
    galleryViewToggle.className = 'gallery-view-toggle';
    galleryViewToggle.innerHTML = '<i class="fas fa-th-large"></i>';
    galleryViewToggle.title = 'Toggle Gallery View';
    document.querySelector('.gallery h2').appendChild(galleryViewToggle);
    
    let isGridView = true;
    galleryViewToggle.addEventListener('click', function() {
        document.querySelector('.gallery-items').classList.toggle('list-view');
        isGridView = !isGridView;
        galleryViewToggle.innerHTML = isGridView ? '<i class="fas fa-th-large"></i>' : '<i class="fas fa-list"></i>';
    });
    
    // Add clear gallery button
    const clearGalleryBtn = document.createElement('button');
    clearGalleryBtn.className = 'clear-gallery-btn';
    clearGalleryBtn.innerHTML = '<i class="fas fa-trash"></i> Clear Gallery';
    clearGalleryBtn.title = 'Clear Gallery';
    document.querySelector('.gallery h2').appendChild(clearGalleryBtn);
    
    clearGalleryBtn.addEventListener('click', function() {
        if (confirm('Are you sure you want to clear your gallery? This cannot be undone.')) {
            galleryItems.innerHTML = '';
        }
    });
    
    // Add tooltip functionality
    const tooltip = document.createElement('div');
    tooltip.className = 'tooltip';
    document.body.appendChild(tooltip);
    
    document.querySelectorAll('[title]').forEach(element => {
        element.addEventListener('mouseover', function(e) {
            tooltip.textContent = this.title;
            tooltip.style.display = 'block';
            tooltip.style.left = e.pageX + 10 + 'px';
            tooltip.style.top = e.pageY + 10 + 'px';
        });
        
        element.addEventListener('mouseout', function() {
            tooltip.style.display = 'none';
        });
        
        element.addEventListener('mousemove', function(e) {
            tooltip.style.left = e.pageX + 10 + 'px';
            tooltip.style.top = e.pageY + 10 + 'px';
        });
    });
    
    // Add welcome message that fades out
    const welcomeMessage = document.createElement('div');
    welcomeMessage.className = 'welcome-message';
    welcomeMessage.innerHTML = `
        <h2>Welcome to Generative Art Gallery</h2>
        <p>Create beautiful algorithmic art with just a few clicks.</p>
        <p>Use the controls above to customize your artwork.</p>
        <button class="btn close-welcome">Get Started</button>
    `;
    document.body.appendChild(welcomeMessage);
    
    document.querySelector('.close-welcome').addEventListener('click', function() {
        welcomeMessage.classList.add('fade-out');
        setTimeout(() => {
            document.body.removeChild(welcomeMessage);
        }, 500);
    });
    
    // Auto-hide welcome message after 5 seconds
    setTimeout(() => {
        if (document.body.contains(welcomeMessage)) {
            welcomeMessage.classList.add('fade-out');
            setTimeout(() => {
                if (document.body.contains(welcomeMessage)) {
                    document.body.removeChild(welcomeMessage);
                }
            }, 500);
        }
    }, 5000);
    
    // Debug info to help troubleshoot blank screen issues
    console.log("DOM fully loaded");
    console.log("Attempting to generate initial artwork...");
});
