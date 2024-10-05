import React, { Component } from 'react';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom';

import Menu from './components/Menu';
import Home from './pages/Home/index';
import RentBuy from './pages/RentBuy';
import About from './pages/About';
import Notfound from './pages/NotFound';

import './css/App.css';


class App extends Component {
  render() {
    return (
      <React.Fragment>
        <Router>
          <div>
            <Switch>
              <Route exact path="/">
                <Home />
              </Route>
              <Route path="/rent-or-buy">
                <Menu />
                <RentBuy />
              </Route>
              <Route path="/about">
                <Menu />
                <About />
              </Route>
              <Route>
                <Notfound />
              </Route>
            </Switch>
          </div>
        </Router>
      </React.Fragment>
    );
  }
}

export default App;
