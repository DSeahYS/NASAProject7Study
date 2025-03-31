// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', function() {
    // Create twinkling stars
    createStars();
    
    // Load recent alerts
    loadRecentAlerts();
    
    // Setup form submission
    document.getElementById('quick-query-form')
        .addEventListener('submit', handleQuickQuery);
});

// Create star background
function createStars() {
    const spaceBg = document.getElementById('space-bg');
    for (let i = 0; i < 100; i++) {
        const star = document.createElement('div');
        star.className = 'star';
        star.style.width = `${Math.random() * 3 + 1}px`;
        star.style.height = star.style.width;
        star.style.left = `${Math.random() * 100}%`;
        star.style.top = `${Math.random() * 100}%`;
        star.style.animationDelay = `${Math.random() * 5}s`;
        star.style.opacity = Math.random() * 0.5 + 0.1;
        spaceBg.appendChild(star);
    }
}

// Load recent alerts from backend
async function loadRecentAlerts() {
    const alertList = document.getElementById('alert-list');
    
    try {
        // In a real implementation, this would fetch from your API
        // const response = await fetch('/api/alerts/recent');
        // const alerts = await response.json();
        
        // Mock data for demonstration
        const alerts = [
            {
                id: 'GRB240331A',
                type: 'GRB',
                time: new Date().toISOString(),
                ra: '08h02m25.44s',
                dec: '+40d51m25.5s',
                significance: 'High'
            },
            {
                id: 'GW240330B',
                type: 'GW',
                time: new Date(Date.now() - 3600000).toISOString(),
                ra: '12h30m45.67s',
                dec: '-05d12m34.5s',
                significance: 'Medium'
            },
            {
                id: 'GRB240330C',
                type: 'GRB',
                time: new Date(Date.now() - 86400000).toISOString(),
                ra: '23h15m12.34s',
                dec: '+15d45m32.1s',
                significance: 'Low'
            }
        ];
        
        // Clear loading message
        alertList.innerHTML = '';
        
        // Add alerts to the list
        alerts.forEach(alert => {
            const alertItem = document.createElement('div');
            alertItem.className = 'alert-item';
            alertItem.innerHTML = `
                <h4>${alert.type} ${alert.id}</h4>
                <p>Time: ${new Date(alert.time).toLocaleString()}</p>
                <p>Coordinates: ${alert.ra}, ${alert.dec}</p>
                <p>Significance: ${alert.significance}</p>
            `;
            alertList.appendChild(alertItem);
        });
        
    } catch (error) {
        alertList.innerHTML = `
            <div class="alert-item error">
                <h4>Error loading alerts</h4>
                <p>${error.message}</p>
            </div>
        `;
        console.error('Error loading alerts:', error);
    }
}

// Handle quick query form submission
async function handleQuickQuery(event) {
    event.preventDefault();
    
    const objectId = document.getElementById('object-id').value;
    const radius = document.getElementById('search-radius').value;
    
    if (!objectId) {
        alert('Please enter an object ID or coordinates');
        return;
    }
    
    try {
        // In a real implementation, this would submit to your API
        // const response = await fetch('/api/query', {
        //     method: 'POST',
        //     headers: { 'Content-Type': 'application/json' },
        //     body: JSON.stringify({ objectId, radius })
        // });
        // const results = await response.json();
        
        // For now, just show a message
        alert(`Query submitted for ${objectId} with radius ${radius}°`);
        
        // In a real app, you would redirect to query results page
        // window.location.href = `query.html?object=${encodeURIComponent(objectId)}&radius=${radius}`;
        
    } catch (error) {
        alert(`Query failed: ${error.message}`);
        console.error('Query error:', error);
    }
}

// Update system status (mock for now)
function updateSystemStatus() {
    // In a real implementation, this would poll your API
    // for system status updates
    setInterval(() => {
        // Randomly toggle status for demo
        const statusItems = document.querySelectorAll('.status-item');
        statusItems.forEach(item => {
            if (Math.random() > 0.8) {
                item.classList.toggle('warning');
                item.classList.toggle('online');
            }
        });
    }, 5000);
}

// Initialize status updates
updateSystemStatus();