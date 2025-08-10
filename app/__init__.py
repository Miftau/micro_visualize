from flask import Flask, Blueprint
from flask_wtf.csrf import CSRFProtect
import os
from wtforms.widgets import html_params
from markupsafe import Markup


def render_field(field, **kwargs):
    return Markup(f"""
        <div class="form-group">
            {field.label}
            {field(**kwargs)}
            {''.join(f'<div class="text-danger">{e}</div>' for e in field.errors)}
        </div>
    """)


def create_app():
    app = Flask(__name__)
    app.jinja_env.globals['render_field'] = render_field
    app.config['SECRET_KEY'] = os.environ.get("SECRET_KEY", "dev_secret_key")
    app.config['UPLOAD_FOLDER'] = os.path.join(os.path.dirname(__file__), 'Uploads')

    # Ensure upload folder exists
    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

    app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024  # 50MB limit

    csrf = CSRFProtect()
    csrf.init_app(app)

    from .routes import bp
    app.register_blueprint(bp)

    return app
