# -*- coding: utf-8 -*-

# Copyright © 2018-2025 Dynare Team
#
# This file is part of Dynare.
#
# Dynare is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Dynare is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
from __future__ import unicode_literals

import os
import sys

sys.path.insert(0, os.path.abspath('../utils'))

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.mathjax',
              'sphinxcontrib.bibtex']

bibtex_bibfiles = ['../../dr.bib']

# define author_year_round style, see https://github.com/mcmtroffaes/sphinxcontrib-bibtex/blob/develop/test/roots/test-citation_style_round_brackets/conf.py
import dataclasses

from dataclasses import dataclass, field
import sphinxcontrib.bibtex.plugin

from sphinxcontrib.bibtex.style.referencing import BracketStyle
from sphinxcontrib.bibtex.style.referencing.author_year import AuthorYearReferenceStyle

def bracket_style() -> BracketStyle:
    return BracketStyle(
        left="(",
        right=")",
    )

@dataclass
class MyReferenceStyle(AuthorYearReferenceStyle):
    bracket_parenthetical: BracketStyle = field(default_factory=bracket_style)
    bracket_textual: BracketStyle = field(default_factory=bracket_style)
    bracket_author: BracketStyle = field(default_factory=bracket_style)
    bracket_label: BracketStyle = field(default_factory=bracket_style)
    bracket_year: BracketStyle = field(default_factory=bracket_style)

sphinxcontrib.bibtex.plugin.register_plugin(
    "sphinxcontrib.bibtex.style.referencing", "author_year_round", MyReferenceStyle
)

# define mystyle for bibliography formatting, based on unsrt.bst from pybtex

# Copyright (c) 2006-2021  Andrey Golovizin

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import re
import pybtex.plugin

from pybtex.richtext import Symbol, Text
from pybtex.style.formatting import BaseStyle, toplevel
from pybtex.style.template import (
    field, first_of, href, join, names, optional, optional_field, sentence,
    tag, together, words
)

def dashify(text):
    dash_re = re.compile(r'-+')
    return Text(Symbol('ndash')).join(text.split(dash_re))


pages = field('pages', apply_func=dashify)

# Keep the year stripped so parentheses render as "(1998)" instead of "( 1998 )".
# Some backends pass a rich Text object; coerce to str before stripping to avoid attribute errors.
date = field('year', apply_func=lambda y: str(y).strip())
# Prebuild a no-space parenthesized year to avoid template spacing adding gaps.
date_in_parens = field('year', apply_func=lambda y: f"({str(y).strip()})")

class MyStyle(BaseStyle):
    default_sorting_style = 'author_year_title'
    
    def format_names(self, role, as_sentence=True):
        formatted_names = names(role, sep=', ', sep2 = ' and ', last_sep=', and ')
        if as_sentence:
            return sentence [formatted_names , date_in_parens]
        else:
            return formatted_names

    def get_article_template(self, e):
        volume_and_pages = first_of [
            # volume and pages, with optional issue number
            optional [
                join [
                    field('volume'),
                    optional['(', field('number'),')'],
                    ':', pages
                ],
            ],
            # pages only
            words ['pages', pages],
        ]
        template = toplevel [
            self.format_names('author'),
            self.format_title(e, 'title'),
            sentence [
                tag('em') [field('journal')],
                optional[ volume_and_pages ]],
            sentence [ optional_field('note') ],
            self.format_web_refs(e),
        ]
        return template

    def format_author_or_editor(self, e):
        return first_of [
            optional[ self.format_names('author') ],
            self.format_editor(e),
        ]

    def format_editor(self, e, as_sentence=True):
        editors = self.format_names('editor', as_sentence=False)
        if 'editor' not in e.persons:
            # when parsing the template, a FieldIsMissing exception
            # will be thrown anyway; no need to do anything now,
            # just return the template that will throw the exception
            return editors
        if len(e.persons['editor']) > 1:
            word = 'editors'
        else:
            word = 'editor'
        result = join(sep=', ') [editors, word]
        if as_sentence:
            return sentence [result]
        else:
            return result

    def format_volume_and_series(self, e, as_sentence=True):
        volume_and_series = optional [
            words [
                together ['Volume' if as_sentence else 'volume', field('volume')], optional [
                    words ['of', field('series')]
                ]
            ]
        ]
        number_and_series = optional [
            words [
                join(sep=Symbol('nbsp')) ['Number' if as_sentence else 'number', field('number')],
                optional [
                    words ['in', field('series')]
                ]
            ]
        ]
        series = optional_field('series')
        result = first_of [
            volume_and_series,
            number_and_series,
            series,
        ]
        if as_sentence:
            return sentence(capfirst=True) [result]
        else:
            return result

    def format_chapter_and_pages(self, e):
        return join(sep=', ') [
            optional [together ['chapter', field('chapter')]],
            optional [together ['pages', pages]],
        ]

    def format_edition(self, e):
        return optional [
            words [
                field('edition', apply_func=lambda x: x.lower()),
                'edition',
            ]
        ]

    def format_title(self, e, which_field, as_sentence=True):
        formatted_title = field(
            which_field, apply_func=lambda text: text.capitalize()
        )
        if as_sentence:
            return sentence [ formatted_title ]
        else:
            return formatted_title

    def format_btitle(self, e, which_field, as_sentence=True):
        formatted_title = tag('em') [ field(which_field) ]
        if as_sentence:
            return sentence[ formatted_title ]
        else:
            return formatted_title

    def format_address_organization_publisher_date(
        self, e, include_organization=True):
        """Format address, organization, publisher, and date.
        Everything is optional, except the date.
        """
        # small difference from unsrt.bst here: unsrt.bst
        # starts a new sentence only if the address is missing;
        # for simplicity here we always start a new sentence
        if include_organization:
            organization = optional_field('organization')
        else:
            organization = None
        return first_of[
            # this will be rendered if there is an address
            optional [
                join(sep=' ') [
                    sentence[
                        field('address'),
                        date,
                    ],
                    sentence[
                        organization,
                        optional_field('publisher'),
                    ],
                ],
            ],
            # if there is no address then we have this
            sentence[
                organization,
                optional_field('publisher'),
                date,
            ],
        ]

    def get_book_template(self, e):
        template = toplevel [
            self.format_author_or_editor(e),
            self.format_btitle(e, 'title'),
            self.format_volume_and_series(e),
            sentence [
                field('publisher'),
                optional_field('address'),
                self.format_edition(e)
            ],
            optional[ sentence [ self.format_isbn(e) ] ],
            sentence [ optional_field('note') ],
            self.format_web_refs(e),
        ]
        return template

    def get_booklet_template(self, e):
        template = toplevel [
            self.format_names('author'),
            self.format_title(e, 'title'),
            sentence [
                optional_field('howpublished'),
                optional_field('address'),
                date,
                optional_field('note'),
            ],
            self.format_web_refs(e),
        ]
        return template

    def get_dataset_template(self, e):
        return self.get_misc_template(e)

    def get_inbook_template(self, e):
        template = toplevel [
            self.format_author_or_editor(e),
            sentence [
                self.format_btitle(e, 'title', as_sentence=False),
                self.format_chapter_and_pages(e),
            ],
            self.format_volume_and_series(e),
            sentence [
                field('publisher'),
                optional_field('address'),
                optional [
                    words [field('edition'), 'edition']
                ],
                date,
                optional_field('note'),
            ],
            self.format_web_refs(e),
        ]
        return template

    def get_incollection_template(self, e):
        template = toplevel [
            sentence [self.format_names('author')],
            self.format_title(e, 'title'),
            words [
                'In',
                sentence [
                    optional[ self.format_editor(e, as_sentence=False) ],
                    self.format_btitle(e, 'booktitle', as_sentence=False),
                    self.format_volume_and_series(e, as_sentence=False),
                    self.format_chapter_and_pages(e),
                ],
            ],
            sentence [
                optional_field('publisher'),
                optional_field('address'),
                self.format_edition(e)
            ],
            self.format_web_refs(e),
        ]
        return template

    def get_inproceedings_template(self, e):
        template = toplevel [
            sentence [self.format_names('author')],
            self.format_title(e, 'title'),
            words [
                'In',
                sentence [
                    optional[ self.format_editor(e, as_sentence=False) ],
                    self.format_btitle(e, 'booktitle', as_sentence=False),
                    self.format_volume_and_series(e, as_sentence=False),
                    optional[ pages ],
                ],
                self.format_address_organization_publisher_date(e),
            ],
            sentence [ optional_field('note') ],
            self.format_web_refs(e),
        ]
        return template

    def get_manual_template(self, e):
        # TODO this only corresponds to the bst style if author is non-empty
        # for empty author we should put the organization first
        template = toplevel [
            optional [ sentence [ self.format_names('author') ] ],
            self.format_btitle(e, 'title'),
            sentence [
                optional_field('organization'),
                optional_field('address'),
                self.format_edition(e),
                optional[ date ],
            ],
            sentence [ optional_field('note') ],
            self.format_web_refs(e),
        ]
        return template

    def get_mastersthesis_template(self, e):
        template = toplevel [
            sentence [self.format_names('author')],
            self.format_title(e, 'title'),
            sentence[
                "Master's thesis",
                field('school'),
                optional_field('address'),
                date,
            ],
            sentence [ optional_field('note') ],
            self.format_web_refs(e),
        ]
        return template

    def get_misc_template(self, e):
        template = toplevel [
            optional[ sentence [self.format_names('author')] ],
            optional[ self.format_title(e, 'title') ],
            sentence[
                optional[ field('howpublished') ],
                optional[ date ],
            ],
            sentence [ optional_field('note') ],
            self.format_web_refs(e),
        ]
        return template

    def get_online_template(self, e):
        return self.get_misc_template(e)

    def get_patent_template(self, e):
        return self.get_misc_template(e)

    def get_phdthesis_template(self, e):
        template = toplevel [
            sentence [self.format_names('author')],
            self.format_btitle(e, 'title'),
            sentence[
                first_of [
                    optional_field('type'),
                    'PhD thesis',
                ],
                field('school'),
                optional_field('address'),
                date,
            ],
            sentence [ optional_field('note') ],
            self.format_web_refs(e),
        ]
        return template

    def get_proceedings_template(self, e):
        if 'editor' in e.persons:
            main_part = [
                self.format_editor(e),
                sentence [
                    self.format_btitle(e, 'title', as_sentence=False),
                    self.format_volume_and_series(e, as_sentence=False),
                    self.format_address_organization_publisher_date(e),
                ],
            ]
        else:
            main_part = [
                optional [ sentence [ field('organization') ] ],
                sentence [
                    self.format_btitle(e, 'title', as_sentence=False),
                    self.format_volume_and_series(e, as_sentence=False),
                    self.format_address_organization_publisher_date(
                        e, include_organization=False),
                ],
            ]
        template = toplevel [
            main_part + [
                sentence [ optional_field('note') ],
                self.format_web_refs(e),
            ]
        ]
        return template

    def get_software_template(self, e):
        return self.get_misc_template(e)

    def get_techreport_template(self, e):
        template = toplevel [
            sentence [self.format_names('author')],
            self.format_title(e, 'title'),
            sentence [
                words[
                    first_of [
                        optional_field('type'),
                        'Technical Report',
                    ],
                    optional_field('number'),
                ],
                field('institution'),
                optional_field('address')
            ],
            sentence [ optional_field('note') ],
            self.format_web_refs(e),
        ]
        return template

    def get_unpublished_template(self, e):
        template = toplevel [
            sentence [self.format_names('author')],
            self.format_title(e, 'title'),
            sentence [
                field('note')
            ],
            self.format_web_refs(e),
        ]
        return template

    def format_web_refs(self, e):
        # based on urlbst output.web.refs
        return sentence [
            optional [ self.format_url(e),
                       optional [ ' (visited on ', field('urldate'), ')' ] ],
            optional [ self.format_eprint(e) ],
            optional [ self.format_pubmed(e) ],
            optional [ self.format_doi(e) ],
            ]

    def format_url(self, e):
        # based on urlbst format.url
        return words [
            'URL:',
            href [
                field('url', raw=True),		
                field('url', raw=True)		
                ]
        ]

    def format_pubmed(self, e):
        # based on urlbst format.pubmed
        return href [            
            join [
                'https://www.ncbi.nlm.nih.gov/pubmed/',		
                field('pubmed', raw=True)		
                ],
            join [
                'PMID:',
                field('pubmed', raw=True)
                ]
            ]

    def format_doi(self, e):
        # based on urlbst format.doi
        return href [
            join [
                'https://doi.org/',		
                field('doi', raw=True)		
                ],
            join [
                'doi:',
                field('doi', raw=True)
                ]
            ]

    def format_eprint(self, e):
        # based on urlbst format.eprint
        return href [
            join [
                'https://arxiv.org/abs/',		
                field('eprint', raw=True)		
                ],
            join [
                'arXiv:',
                field('eprint', raw=True)
                ]
            ]

    def format_isbn(self, e):
        return join(sep=' ') [ 'ISBN', field('isbn') ]


#end definition mystyle 

pybtex.plugin.register_plugin("pybtex.style.formatting", "mystyle", MyStyle)


bibtex_default_style = 'mystyle'

bibtex_reference_style = 'author_year_round'

source_suffix = '.rst'

templates_path = ['_templates']

html_static_path = ['_static']

master_doc = 'index'

project = u'Dynare'
copyright = u'1996–2025 Dynare Team'
author = u'Dynare Team'

add_function_parentheses = False

language = 'en'

exclude_patterns = []

highlight_language = 'dynare'

todo_include_todos = False

html_theme = 'alabaster'

html_sidebars = {
    "**": [
        "about.html",
        "searchbox.html",
        "navigation.html",
    ]
}

html_theme_options = {
    'logo': 'dlogo.svg',
    'logo_name': False,
    'fixed_sidebar': True,
    'page_width': '100%',
}

htmlhelp_basename = 'Dynaremanual'

latex_elements = {
    'sphinxsetup': 'VerbatimBorderColor={rgb}{1,1,1},VerbatimColor={RGB}{240,240,240}, \
                    warningBorderColor={RGB}{255,50,50},OuterLinkColor={RGB}{34,139,34}, \
                    InnerLinkColor={RGB}{51,51,255},TitleColor={RGB}{51,51,255}',
    'papersize': 'a4paper',
    # Add support for the perpendicular symbol input as UTF-8
    'preamble': r'''
\DeclareUnicodeCharacter{27C2}{\ensuremath{\perp}}
'''
}

latex_documents = [
    (master_doc, 'dynare-manual.tex', u'Dynare Reference Manual',
     u'Dynare Team', 'manual'),
]

man_pages = [
    (master_doc, 'dynare', u'Dynare Reference Manual',
     [author], 1)
]

def setup(app):
    from dynare_dom import DynareDomain
    from dynare_lex import DynareLexer
    app.add_lexer("dynare", DynareLexer)
    app.add_domain(DynareDomain)
