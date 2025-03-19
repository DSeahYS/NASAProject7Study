from flask import Flask, render_template, jsonify, request, send_file
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Set the backend to Agg (important for headless environments)
import matplotlib.pyplot as plt
import random
import math
import os
import io
import base64
from matplotlib.colors import LinearSegmentedColormap
from datetime import datetime

app = Flask(__name__)

class GenerativeArtCreator:
    def __init__(self):
        self.data_x = None
        self.data_y = None
        self.data_z = None
        self.seed = None
    
    def generate_art(self, style="spiral", seed=None, points=50000):
        """Generate art based on the selected style"""
        if style == "spiral":
            return self.generate_spiral_waves(seed, points)
        elif style == "fractal":
            return self.generate_fractal_dream(seed, points)
        elif style == "cosmic":
            return self.generate_cosmic_flow(seed, points)
        elif style == "harmonic":
            return self.generate_harmonic_structure(seed, points)
        else:
            return self.generate_spiral_waves(seed, points)
    
    def generate_spiral_waves(self, seed=None, points=50000):
        """Generate spiral wave pattern"""
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)
        else:
            self.seed = random.randint(1, 10000)
            random.seed(self.seed)
            np.random.seed(self.seed)
        
        self.data_x = []
        self.data_y = []
        self.data_z = []
        
        for _ in range(points):
            a = random.uniform(0, 2*math.pi)
            r = random.uniform(0, 1)
            
            x = r * math.cos(a + r * 10) * math.sin(r * 5)
            y = r * math.sin(a + r * 10) * math.sin(r * 5)
            z = random.random()
            
            self.data_x.append(x)
            self.data_y.append(y)
            self.data_z.append(z)
        
        return self
    
    def generate_fractal_dream(self, seed=None, points=50000):
        """Generate fractal-like pattern"""
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)
        else:
            self.seed = random.randint(1, 10000)
            random.seed(self.seed)
            np.random.seed(self.seed)
        
        self.data_x = []
        self.data_y = []
        self.data_z = []
        
        for _ in range(points):
            x = random.uniform(-1, 1)
            y = random.uniform(-1, 1)
            
            for _ in range(5):
                nx = math.sin(y * 3) - 0.2 * math.sin(x * 2)
                ny = math.sin(x * 3) - 0.2 * math.sin(y * 2)
                x, y = nx, ny
            
            z = (x**2 + y**2) ** 0.5
            
            self.data_x.append(x)
            self.data_y.append(y)
            self.data_z.append(z)
        
        return self
    
    def generate_cosmic_flow(self, seed=None, points=50000):
        """Generate flowing cosmic-like pattern"""
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)
        else:
            self.seed = random.randint(1, 10000)
            random.seed(self.seed)
            np.random.seed(self.seed)
        
        self.data_x = []
        self.data_y = []
        self.data_z = []
        
        for _ in range(points):
            t = random.uniform(0, 20 * math.pi)
            scale = random.uniform(0.1, 1.0)
            
            x = math.sin(t) * (math.exp(math.cos(t)) - 2 * math.cos(4*t) - math.sin(t/12)**5) * scale
            y = math.cos(t) * (math.exp(math.cos(t)) - 2 * math.cos(4*t) - math.sin(t/12)**5) * scale
            z = t / (20 * math.pi)
            
            self.data_x.append(x)
            self.data_y.append(y)
            self.data_z.append(z)
        
        return self
    
    def generate_harmonic_structure(self, seed=None, points=50000):
        """Generate harmonic structure pattern"""
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)
        else:
            self.seed = random.randint(1, 10000)
            random.seed(self.seed)
            np.random.seed(self.seed)
        
        self.data_x = []
        self.data_y = []
        self.data_z = []
        
        for _ in range(points):
            theta = random.uniform(0, 2 * math.pi)
            phi = random.uniform(0, 2 * math.pi)
            r = random.uniform(0, 1)
            
            x = r * math.sin(theta) * math.cos(phi + r * 5)
            y = r * math.sin(theta) * math.sin(phi + r * 5)
            z = math.cos(theta + r * 3)
            
            self.data_x.append(x)
            self.data_y.append(y)
            self.data_z.append(z)
        
        return self
    
    def get_image_data(self, style="dark", cmap_name="cosmic", size=10, dpi=100, point_size=0.5, alpha=0.7):
        """Generate image data as base64 string"""
        if self.data_x is None or self.data_y is None:
            return None
        
        # Define color maps
        color_maps = {
            "cosmic": LinearSegmentedColormap.from_list("cosmic", 
                     [(0, 0, 0.3), (0.1, 0, 0.6), (0.3, 0, 0.9), 
                      (0.5, 0.2, 1), (0.7, 0.6, 1), (1, 1, 1)]),
            "sunset": LinearSegmentedColormap.from_list("sunset", 
                     [(0.1, 0, 0.2), (0.3, 0, 0.5), (0.5, 0.2, 0.4), 
                      (0.7, 0.6, 0.2), (0.9, 0.9, 0.1), (1, 1, 0.8)]),
            "ocean": LinearSegmentedColormap.from_list("ocean", 
                    [(0, 0, 0.3), (0, 0.3, 0.6), (0, 0.6, 0.8), 
                     (0.1, 0.8, 0.9), (0.5, 1, 1)]),
            "forest": LinearSegmentedColormap.from_list("forest", 
                     [(0, 0.1, 0), (0.1, 0.3, 0.1), (0.3, 0.5, 0.2), 
                      (0.5, 0.7, 0.3), (0.7, 0.9, 0.5), (0.9, 1, 0.7)])
        }
        
        plt.figure(figsize=(size, size), dpi=dpi, facecolor='black' if style == "dark" else 'white')
        
        if cmap_name in color_maps:
            colormap = color_maps[cmap_name]
        else:
            colormap = plt.cm.viridis
        
        plt.scatter(self.data_x, self.data_y, c=self.data_z, cmap=colormap, 
                   s=point_size, alpha=alpha, edgecolors='none')
        
        plt.axis('off')
        if style == "dark":
            plt.gca().set_facecolor('black')
        
        plt.tight_layout(pad=0)
        
        # Convert plot to base64 string
        buffer = io.BytesIO()
        plt.savefig(buffer, format='png', bbox_inches='tight', pad_inches=0)
        buffer.seek(0)
        image_png = buffer.getvalue()
        buffer.close()
        plt.close()
        
        image_data = base64.b64encode(image_png).decode('utf-8')
        return image_data


# Initialize the art creator
art_creator = GenerativeArtCreator()

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/generate', methods=['GET'])
def generate_art():
    style = request.args.get('style', 'spiral')
    cmap = request.args.get('cmap', 'cosmic')
    seed_param = request.args.get('seed')
    
    seed = None
    if seed_param:
        try:
            seed = int(seed_param)
        except ValueError:
            pass
    
    # Generate the art
    art_creator.generate_art(style=style, seed=seed)
    image_data = art_creator.get_image_data(cmap_name=cmap)
    
    return jsonify({
        'image': image_data,
        'seed': art_creator.seed
    })

@app.route('/animate', methods=['GET'])
def animate_art():
    style = request.args.get('style', 'spiral')
    cmap = request.args.get('cmap', 'cosmic')
    frames_param = request.args.get('frames', '30')
    
    try:
        frames = int(frames_param)
    except ValueError:
        frames = 30
    
    # Create animation frames
    frames_data = []
    base_seed = random.randint(1, 10000)
    
    for i in range(frames):
        # Generate slightly modified art for each frame
        seed = base_seed + i
        art_creator.generate_art(style=style, seed=seed)
        image_data = art_creator.get_image_data(cmap_name=cmap)
        frames_data.append(image_data)
    
    return jsonify({
        'frames': frames_data,
        'base_seed': base_seed
    })

@app.route('/download', methods=['POST'])
def download_art():
    # Get image data from request
    data = request.json
    image_data = data.get('image')
    
    if not image_data:
        return jsonify({'error': 'No image data provided'}), 400
    
    # Convert base64 to image
    if ',' in image_data:
        image_data = image_data.split(',')[1]
    
    try:
        image_bytes = base64.b64decode(image_data)
    except Exception as e:
        return jsonify({'error': f'Invalid image data: {str(e)}'}), 400
    
    # Create a temporary file
    temp_file = io.BytesIO(image_bytes)
    temp_file.seek(0)
    
    # Generate filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"generative_art_{timestamp}.png"
    
    # Send the file
    return send_file(
        temp_file,
        mimetype='image/png',
        as_attachment=True,
        download_name=filename
    )

@app.errorhandler(404)
def page_not_found(e):
    return "Page not found", 404

@app.errorhandler(500)
def internal_server_error(e):
    return "Internal server error", 500

if __name__ == '__main__':
    # Make sure the templates and static folders exist
    if not os.path.exists('templates'):
        os.makedirs('templates')
    if not os.path.exists('static/css'):
        os.makedirs('static/css')
    if not os.path.exists('static/js'):
        os.makedirs('static/js')
    
    app.run(debug=True)
