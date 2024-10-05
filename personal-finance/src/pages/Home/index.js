import React, { Component } from 'react';
import ReactTooltip from 'react-tooltip';

import CashFlow from './CashFlow';
import IncomeExpenses from './IncomeExpenses';
import AssetsLiabilities from './AssetsLiabilities';
import NetWorth from './NetWorth';
import Retirement from './Retirement';

import Header from '../../components/Header';
import Footer from '../../components/Footer';
import Donut from '../../components/Donut';
import Storage from '../../components/Storage';

import CONFIG from '../../config';

class Home extends Component {
  constructor(props) {
    super(props);

    this.state = {
      all_items: Storage.getItems(),
      assets: Storage.getItems(CONFIG.ITEM_TYPE.ASSET),
      libilities: Storage.getItems(CONFIG.ITEM_TYPE.LIABILITIES),
      income: Storage.getItems(CONFIG.ITEM_TYPE.INCOME),
      expenses: Storage.getItems(CONFIG.ITEM_TYPE.EXPENSE),
      sample: false,
      active_tab: 'incomeexpenses',
    };

    this.loadSampleData = this.loadSampleData.bind(this);
  }

  loadSampleData() {
    let all_items = [];

    all_items = all_items.concat(
      CONFIG.SAMPLE_DATA.ASSETS,
      CONFIG.SAMPLE_DATA.LIABILITIES,
      CONFIG.SAMPLE_DATA.INCOME,
      CONFIG.SAMPLE_DATA.EXPENSES
    );

    this.setState({
      sample: true,
      all_items: all_items,
      assets: CONFIG.SAMPLE_DATA.ASSETS,
      libilities: CONFIG.SAMPLE_DATA.LIABILITIES,
      income: CONFIG.SAMPLE_DATA.INCOME,
      expenses: CONFIG.SAMPLE_DATA.EXPENSES,
    });
  }

  changeTab = (tab) => {
    this.setState({
      active_tab: tab,
    });
  };

  renderHeader() {
    return (
      <Header items={this.state.all_items} sampleHandler={this.loadSampleData}>
        {this.state.all_items.length > 0 && (
          <Donut items={this.state.all_items} />
        )}
      </Header>
    );
  }

  renderIncomeExpenses() {
    return (
      <React.Fragment>
        {!this.state.sample && <IncomeExpenses items={this.state.all_items} />}
        {this.state.sample && <IncomeExpenses items={this.state.all_items} />}
      </React.Fragment>
    );
  }

  renderAssetsLiabilities() {
    return (
      <React.Fragment>
        {!this.state.sample && (
          <AssetsLiabilities items={this.state.all_items} />
        )}
        {this.state.sample && (
          <AssetsLiabilities items={this.state.all_items} />
        )}
      </React.Fragment>
    );
  }

  renderCashflow() {
    return (
      <React.Fragment>
        {!this.state.sample && <CashFlow items={this.state.all_items} />}
        {this.state.sample && <CashFlow items={this.state.all_items} />}
      </React.Fragment>
    );
  }

  renderNetWorth() {
    return (
      <React.Fragment>
        {!this.state.sample && <NetWorth items={this.state.all_items} />}
        {this.state.sample && <NetWorth items={this.state.all_items} />}
      </React.Fragment>
    );
  }

  renderRetirement() {
    return (
      <React.Fragment>
        {!this.state.sample && <Retirement items={this.state.all_items} />}
        {this.state.sample && <Retirement items={this.state.all_items} />}
      </React.Fragment>
    );
  }

  render() {
    return (
      <React.Fragment>
        <ReactTooltip effect="solid" />

        {!this.state.sample && this.renderHeader()}
        {this.state.sample && this.renderHeader()}

        <div className="bg-white" style={{ marginTop: '-50px' }}>
          <div className="container">
            <div className="row">
              <div className="col">
                <ul className="nav nav-tabs" role="tablist">
                  <li className="nav-item">
                    <a
                      className="nav-link active"
                      data-toggle="tab"
                      href="javascript:;"
                      onClick={this.changeTab.bind(this, 'incomeexpenses')}
                    >
                      Income &amp; Expenses
                    </a>
                  </li>
                  <li className="nav-item">
                    <a
                      className="nav-link"
                      data-toggle="tab"
                      href="javascript:;"
                      onClick={this.changeTab.bind(this, 'assetsliabilities')}
                    >
                      Assets &amp; Liabilities
                    </a>
                  </li>
                  <li className="nav-item">
                    <a
                      className="nav-link"
                      data-toggle="tab"
                      href="javascript:;"
                      onClick={this.changeTab.bind(this, 'cashflow')}
                    >
                      Cash Flow
                    </a>
                  </li>
                  <li className="nav-item">
                    <a
                      className="nav-link"
                      data-toggle="tab"
                      href="javascript:;"
                      onClick={this.changeTab.bind(this, 'networth')}
                    >
                      Net Worth
                    </a>
                  </li>
                  <li className="nav-item">
                    <a
                      className="nav-link"
                      data-toggle="tab"
                      href="javascript:;"
                      onClick={this.changeTab.bind(this, 'retirement')}
                    >
                      Retirement
                    </a>
                  </li>
                </ul>
              </div>
            </div>
          </div>
        </div>

        <br />

        {this.state.active_tab === 'incomeexpenses' &&
          this.renderIncomeExpenses()}
        {this.state.active_tab === 'assetsliabilities' &&
          this.renderAssetsLiabilities()}
        {this.state.active_tab === 'cashflow' && this.renderCashflow()}
        {this.state.active_tab === 'networth' && this.renderNetWorth()}
        {this.state.active_tab === 'retirement' && this.renderRetirement()}

        <Footer />
      </React.Fragment>
    );
  }
}

export default Home;
