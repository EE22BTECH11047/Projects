from flask import Flask, send_file, render_template, jsonify
import random
import os
import base64

app = Flask(__name__, template_folder='/storage/self/primary/Download/EE23010/audio-player/src', static_folder='/storage/self/primary/Download/EE23010/audio-player/src')

audio_directory = '/storage/self/primary/Download/EE23010/audio-player/audio'

@app.route('/play-random-audio')
def play_random():
    audio_files = os.listdir(audio_directory)
    return render_template('pages/front.html', audio_files=audio_files)

@app.route('/get-shuffled-indices')
def get_shuffled_indices():
    audio_files = os.listdir(audio_directory)
    uni_dist = [random.uniform(0, 1) for _ in range(len(audio_files))]
    indices = [int(rand_no * len(audio_files)) for rand_no in uni_dist]
    return indices

@app.route('/get-audio-list')
def get_audio_list():
    audio_files = os.listdir(audio_directory)
    return jsonify(audio_files)

@app.route('/audio/<filename>')
def serve_audio(filename):
    audio_path = os.path.join(audio_directory, filename)
    return send_file(audio_path)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
