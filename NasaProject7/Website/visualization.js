class DataVisualizer {
    constructor() {
        this.plotArea = document.getElementById('plot-area');
        this.plotType = document.getElementById('plot-type');
        this.dataSource = document.getElementById('data-source');
        this.generateBtn = document.getElementById('generate-plot');
        this.showErrorBars = document.getElementById('show-error-bars');
        this.showModel = document.getElementById('show-model');
        this.colorScheme = document.getElementById('color-scheme');
        
        this.currentPlot = null;
        this.setupEventListeners();
    }
    
    setupEventListeners() {
        this.generateBtn.addEventListener('click', () => this.generatePlot());
    }
    
    async generatePlot() {
        try {
            const plotType = this.plotType.value;
            const dataSource = this.dataSource.value;
            
            // Clear previous plot
            this.plotArea.innerHTML = '<div class="plot-loading">Generating plot...</div>';
            
            // In a real app, this would fetch data from your API
            // const response = await fetch(`/api/visualization?type=${plotType}&source=${dataSource}`);
            // const plotData = await response.json();
            
            // Generate mock data based on selected plot type
            const plotData = this.generateMockData(plotType);
            
            // Create the appropriate plot
            this.createPlot(plotType, plotData);
            
        } catch (error) {
            this.showError('Failed to generate plot', error);
        }
    }
    
    generateMockData(plotType) {
        // Generate mock data based on plot type
        switch(plotType) {
            case 'light-curve':
                return this.generateLightCurveData();
            case 'finding-chart':
                return this.generateFindingChartData();
            case 'color-magnitude':
                return this.generateColorMagData();
            default:
                throw new Error('Unknown plot type');
        }
    }
    
    generateLightCurveData() {
        const data = {
            times: [],
            magnitudes: [],
            errors: []
        };
        
        // Generate 20 data points over 10 days
        for (let i = 0; i < 20; i++) {
            data.times.push(i * 0.5); // Days
            data.magnitudes.push(16 + Math.sin(i/3) * 2 + Math.random() * 0.3);
            data.errors.push(0.1 + Math.random() * 0.1);
        }
        
        return data;
    }
    
    generateFindingChartData() {
        const data = {
            ra: [],
            dec: [],
            magnitudes: []
        };
        
        // Generate 50 random sources
        for (let i = 0; i < 50; i++) {
            data.ra.push(120 + (Math.random() - 0.5) * 10); // RA in degrees
            data.dec.push(40 + (Math.random() - 0.5) * 10); // Dec in degrees
            data.magnitudes.push(15 + Math.random() * 5);
        }
        
        // Add a bright transient in the center
        data.ra.push(120);
        data.dec.push(40);
        data.magnitudes.push(12);
        
        return data;
    }
    
    generateColorMagData() {
        const data = {
            colors: [],
            magnitudes: []
        };
        
        // Generate main sequence stars
        for (let i = 0; i < 100; i++) {
            const color = 0.5 + Math.random() * 2;
            data.colors.push(color);
            data.magnitudes.push(10 + color * 2 + (Math.random() - 0.5));
        }
        
        // Add some outliers
        for (let i = 0; i < 5; i++) {
            data.colors.push(1 + Math.random());
            data.magnitudes.push(12 + Math.random() * 3);
        }
        
        return data;
    }
    
    createPlot(plotType, data) {
        // Clear previous plot
        this.plotArea.innerHTML = '';
        
        switch(plotType) {
            case 'light-curve':
                this.createLightCurve(data);
                break;
            case 'finding-chart':
                this.createFindingChart(data);
                break;
            case 'color-magnitude':
                this.createColorMagDiagram(data);
                break;
        }
    }
    
    createLightCurve(data) {
        const ctx = document.createElement('canvas');
        this.plotArea.appendChild(ctx);
        
        new Chart(ctx, {
            type: 'line',
            data: {
                labels: data.times.map(t => t.toFixed(1) + 'd'),
                datasets: [{
                    label: 'Magnitude',
                    data: data.magnitudes,
                    borderColor: '#00E5FF',
                    backgroundColor: 'rgba(0, 229, 255, 0.1)',
                    borderWidth: 2,
                    pointRadius: 4,
                    pointBackgroundColor: '#00E5FF',
                    pointBorderColor: '#fff',
                    pointHoverRadius: 6,
                    fill: true,
                    tension: 0.1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    y: {
                        reverse: true,
                        title: {
                            display: true,
                            text: 'Magnitude'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Time (days)'
                        }
                    }
                },
                plugins: {
                    legend: {
                        display: false
                    }
                }
            }
        });
    }
    
    createFindingChart(data) {
        const plotData = [{
            type: 'scatter',
            mode: 'markers',
            x: data.ra,
            y: data.dec,
            marker: {
                size: data.magnitudes.map(m => 20 - m), // Bigger for brighter stars
                color: data.magnitudes,
                colorscale: 'Viridis',
                reversescale: true,
                colorbar: {
                    title: 'Magnitude'
                }
            }
        }];
        
        const layout = {
            title: 'Finding Chart',
            xaxis: {
                title: 'Right Ascension (deg)',
                scaleanchor: 'y'
            },
            yaxis: {
                title: 'Declination (deg)'
            },
            height: 500,
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: {
                color: '#E0E0E0'
            }
        };
        
        Plotly.newPlot(this.plotArea, plotData, layout);
    }
    
    createColorMagDiagram(data) {
        const ctx = document.createElement('canvas');
        this.plotArea.appendChild(ctx);
        
        new Chart(ctx, {
            type: 'scatter',
            data: {
                datasets: [{
                    label: 'Stars',
                    data: data.colors.map((c, i) => ({x: c, y: data.magnitudes[i]})),
                    backgroundColor: '#00E5FF',
                    borderColor: '#fff',
                    borderWidth: 1,
                    pointRadius: 4
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    y: {
                        reverse: true,
                        title: {
                            display: true,
                            text: 'Magnitude'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Color Index'
                        }
                    }
                },
                plugins: {
                    legend: {
                        display: false
                    }
                }
            }
        });
    }
    
    showError(message, error) {
        console.error(message, error);
        this.plotArea.innerHTML = `
            <div class="plot-error">
                <h3>${message}</h3>
                <p>${error.message}</p>
            </div>
        `;
    }
}

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    new DataVisualizer();
});