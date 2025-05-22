// Alert Monitor
class AlertMonitor {
    constructor() {
        this.monitoring = false;
        this.alerts = [];
        this.selectedAlert = null;
        this.socket = null;
        
        // Initialize UI elements
        this.alertTable = document.getElementById('alert-table-body');
        this.alertDetails = document.getElementById('alert-details').querySelector('.details-content');
        this.startBtn = document.getElementById('start-monitoring');
        this.pauseBtn = document.getElementById('pause-monitoring');
        this.simulateBtn = document.getElementById('simulate-alerts');
        
        // Set up event listeners
        this.setupEventListeners();
        
        // Load initial alerts
        this.loadInitialAlerts();
    }
    
    setupEventListeners() {
        this.startBtn.addEventListener('click', () => this.startMonitoring());
        this.pauseBtn.addEventListener('click', () => this.pauseMonitoring());
        this.simulateBtn.addEventListener('click', () => this.simulateAlerts());
        
        // Filter checkboxes
        document.querySelectorAll('input[name="alert-type"]').forEach(checkbox => {
            checkbox.addEventListener('change', () => this.filterAlerts());
        });
    }
    
    async loadInitialAlerts() {
            
            this.renderAlerts();
        } catch (error) {
            this.showError('Failed to load alerts', error);
        }
    }
    
    pauseMonitoring() {
        this.monitoring = false;
        this.startBtn.disabled = false;
        this.pauseBtn.disabled = true;
        
        clearInterval(this.alertInterval);
        console.log('Alert monitoring paused');
    }
    
    simulateAlerts() {
        // Generate a few random alerts
        for (let i = 0; i < 3; i++) {
            setTimeout(() => {
                this.receiveAlert(this.generateRandomAlert());
            }, i * 1000);
        }
    }
    
    generateRandomAlert() {
        const types = ['GRB', 'GW', 'Test'];
        const sigs = ['Low', 'Medium', 'High'];
        const now = new Date();
        
        return {
            id: `${types[Math.floor(Math.random() * types.length)]}${now.getFullYear()}${(now.getMonth()+1).toString().padStart(2,'0')}${now.getDate().toString().padStart(2,'0')}${String.fromCharCode(65+Math.floor(Math.random()*26))}`,
            type: types[Math.floor(Math.random() * types.length)],
            time: now.toISOString(),
            ra: `${Math.floor(Math.random()*24).toString().padStart(2,'0')}h${Math.floor(Math.random()*60).toString().padStart(2,'0')}m${Math.floor(Math.random()*60).toString().padStart(2,'0')}s`,
            dec: `${Math.floor(Math.random()*90)}d${Math.floor(Math.random()*60).toString().padStart(2,'0')}m${Math.floor(Math.random()*60).toString().padStart(2,'0')}s`,
            significance: sigs[Math.floor(Math.random() * sigs.length)],
            details: 'Simulated alert for testing purposes'
        };
    }
    
    receiveAlert(alert) {
        this.alerts.unshift(alert); // Add to beginning
        this.renderAlerts();
    }
    
    filterAlerts() {
        const checkedTypes = Array.from(document.querySelectorAll('input[name="alert-type"]:checked'))
            .map(el => el.value);
        
        const filtered = checkedTypes.length === 0 
            ? []
            : this.alerts.filter(alert => checkedTypes.includes(alert.type));
        
        this.renderAlerts(filtered);
    }
    
    renderAlerts(alerts = this.alerts) {
        if (alerts.length === 0) {
            this.alertTable.innerHTML = '<div class="loading">No alerts found</div>';
            return;
        }
        
        this.alertTable.innerHTML = '';
        
        alerts.forEach(alert => {
            const row = document.createElement('div');
            row.className = 'alert-row';
            if (this.selectedAlert && this.selectedAlert.id === alert.id) {
                row.classList.add('selected');
            }
            
            row.innerHTML = `
                <div>${new Date(alert.time).toLocaleString()}</div>
                <div>${alert.type}</div>
                <div>${alert.id}</div>
                <div>${alert.ra}</div>
                <div>${alert.dec}</div>
                <div>${alert.significance}</div>
            `;
            
            row.addEventListener('click', () => this.showAlertDetails(alert));
            this.alertTable.appendChild(row);
        });
    }
    
    showAlertDetails(alert) {
        this.selectedAlert = alert;
        this.alertDetails.innerHTML = `
            <h3>${alert.type} ${alert.id}</h3>
            <p><strong>Time:</strong> ${new Date(alert.time).toLocaleString()}</p>
            <p><strong>Coordinates:</strong> ${alert.ra}, ${alert.dec}</p>
            <p><strong>Significance:</strong> ${alert.significance}</p>
            <p><strong>Details:</strong> ${alert.details}</p>
        `;
        
        // Update selected row styling
        document.querySelectorAll('.alert-row').forEach(row => {
            row.classList.remove('selected');
        });
        document.querySelectorAll('.alert-row').forEach(row => {
            if (row.textContent.includes(alert.id)) {
                row.classList.add('selected');
            }
        });
    }
    
    showError(message, error) {
        console.error(message, error);
        this.alertTable.innerHTML = `
            <div class="alert-item error">
                <h4>${message}</h4>
                <p>${error.message}</p>
            </div>
        `;
    }
}

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    new AlertMonitor();
});