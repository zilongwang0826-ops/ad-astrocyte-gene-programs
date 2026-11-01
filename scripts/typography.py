"""Use installed academic fonts without downloading or redistributing font files."""
from pathlib import Path
import os

FONT_FILES = {
    'Body': 'Times New Roman.ttf',
    'BodyBold': 'Times New Roman Bold.ttf',
    'BodyItalic': 'Times New Roman Italic.ttf',
    'BodyBoldItalic': 'Times New Roman Bold Italic.ttf',
    'Sans': 'Arial.ttf',
    'SansBold': 'Arial Bold.ttf',
    'SansItalic': 'Arial Italic.ttf',
}

def font_directory():
    configured = os.environ.get('AD_ASTRO_FONT_DIR')
    candidates = [Path(configured)] if configured else [Path('/System/Library/Fonts/Supplemental')]
    for directory in candidates:
        if all((directory / name).is_file() for name in FONT_FILES.values()):
            return directory
    raise RuntimeError('Installed Times New Roman and Arial TrueType files are required. '
                       'Set AD_ASTRO_FONT_DIR to a local directory containing the filenames '
                       'listed in scripts/typography.py. No font download is performed.')

def register_report_fonts():
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    directory = font_directory()
    for name, filename in FONT_FILES.items():
        pdfmetrics.registerFont(TTFont(name, str(directory / filename)))
    pdfmetrics.registerFontFamily('Body', normal='Body', bold='BodyBold',
                                 italic='BodyItalic', boldItalic='BodyBoldItalic')
    pdfmetrics.registerFontFamily('Sans', normal='Sans', bold='SansBold',
                                 italic='SansItalic', boldItalic='SansBold')
    return directory

def register_plot_fonts():
    from matplotlib import font_manager
    directory = font_directory()
    for key in ['Sans', 'SansBold', 'SansItalic']:
        font_manager.fontManager.addfont(str(directory / FONT_FILES[key]))
    return directory
