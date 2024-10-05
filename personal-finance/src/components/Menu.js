import React, { Component } from 'react';
import { Link } from 'react-router-dom';

class Menu extends Component {
  render() {
    return (
      <div className="nav-container">
        <div className="bg-dark navbar-dark position-fixed" data-sticky="top">
          <div className="container">
            <nav className="navbar navbar-expand-lg">
              <a className="navbar-brand" href="/">
                <i className="icon-pie-chart" />
              </a>

              <button
                className="navbar-toggler"
                type="button"
                data-toggle="collapse"
                data-target="#navbarNav"
                aria-controls="navbarNav"
                aria-expanded="false"
                aria-label="Toggle navigation"
              >
                <i className="icon-menu h4"></i>
              </button>
              <div
                className="collapse navbar-collapse justify-content-between"
                id="navbarNav"
              >
                <ul className="navbar-nav">
                  <li className="nav-item">
                    <Link to="/" className="nav-link">
                      Dashboard
                    </Link>
                  </li>
                  <li className="nav-item">
                    <Link to="/rent-or-buy" className="nav-link">
                      Rent or Buy
                    </Link>
                  </li>
                  <li className="nav-item">
                    <Link to="/about" className="nav-link">
                      About
                    </Link>
                  </li>
                  {/*
                                    <li className="nav-item dropdown">
                                        <a className="nav-link dropdown-toggle" href="#" id="pagesDropdown" role="button" data-toggle="dropdown">Pages</a>
                                        <div className="dropdown-menu" aria-labelledby="pagesDropdown">
                                            <a className="dropdown-item" href="pages-landing.html">
                                                <span className="h6 mb-0">Landing Pages</span>
                                                <p className="text-small text-muted">Showcase your product in style</p>
                                            </a>

                                            <div className="dropdown-divider"></div>

                                            <a className="dropdown-item" href="pages-app.html">
                                                <span className="h6 mb-0">App Pages</span>
                                                <p className="text-small text-muted">Build detailed, functional web apps</p>
                                            </a>

                                            <div className="dropdown-divider"></div>

                                            <a className="dropdown-item" href="pages-inner.html">
                                                <span className="h6 mb-0">Inner Pages</span>
                                                <p className="text-small text-muted">Complete your online presence</p>
                                            </a>
                                        </div>
                                    </li>
                                    <li className="nav-item dropdown">
                                        <a className="nav-link dropdown-toggle" href="#" id="componentsDropdown" role="button" data-toggle="dropdown">Components</a>
                                        <div className="dropdown-menu" aria-labelledby="componentsDropdown">
                                            <a className="dropdown-item" href="components-bootstrap.html">Bootstrap</a>
                                            <a className="dropdown-item" href="components-wingman.html">Wingman</a>
                                        </div>
                                    </li>
                                    */}
                </ul>

                {/*
                                <ul className="navbar-nav">
                                    <li className="nav-item dropdown">
                                        <a className="nav-link dropdown-toggle dropdown-toggle-no-arrow p-lg-0" href="http://example.com" id="dropdown01" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                                            <img alt="Image" src="assets/img/avatar-male-3.jpg" className="avatar avatar-xs" />
                                            <span className="badge badge-danger">8</span>
                                        </a>
                                        <div className="dropdown-menu dropdown-menu-sm dropdown-menu-right" aria-labelledby="dropdown01">
                                            <a className="dropdown-item" href="#">Notifications <span className="badge badge-danger">8</span></a>
                                            <a className="dropdown-item" href="#">Profile</a>
                                            <div className="dropdown-divider"></div>
                                            <a className="dropdown-item" href="#">Settings</a>
                                            <a className="dropdown-item" href="#">Groups</a>
                                            <a className="dropdown-item" href="#">Log out</a>
                                        </div>
                                    </li>
                                </ul>
                                */}
              </div>
            </nav>
          </div>
        </div>
      </div>
    );
  }
}

export default Menu;
