# microhub_flask/app/__init__.py
from flask import Flask
from flask_wtf.csrf import CSRFProtect
import os

csrf = CSRFProtect()

def create_app():
    app = Flask(__name__)
    app.config['SECRET_KEY'] = os.environ.get("SECRET_KEY", "dev_secret_key")
    app.config['UPLOAD_FOLDER'] = os.path.join(os.getcwd(), 'microhub_flask', 'uploads')
    app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024  # 50MB limit

    csrf.init_app(app)

    from .routes import main
    app.register_blueprint(main)

    return app
