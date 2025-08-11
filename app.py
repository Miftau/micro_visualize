
from flask import Flask, Blueprint
from flask_wtf.csrf import CSRFProtect
import os
from wtforms.widgets import html_params
from markupsafe import Markup


app = Flask(__name__)
app.config['SECRET_KEY'] = os.environ.get("SECRET_KEY", "dev_secret_key")
app.config['UPLOAD_FOLDER'] = os.path.join(os.path.dirname(__file__), 'Uploads')
from server.routes import bp
app.register_blueprint(bp)

def render_field(field, **kwargs):
    return Markup(f"""
        <div class="form-group">
            {field.label}
            {field(**kwargs)}
            {''.join(f'<div class="text-danger">{e}</div>' for e in field.errors)}
        </div>
    """)
app.jinja_env.globals['render_field'] = render_field


# Ensure upload folder exists
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024  # 50MB limit

csrf = CSRFProtect()
csrf.init_app(app)


if __name__ == "__main__":
    app.run()
