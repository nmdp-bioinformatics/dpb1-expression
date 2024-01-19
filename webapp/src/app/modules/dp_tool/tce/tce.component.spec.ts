import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { TceComponent } from './tce.component';

describe('TceComponent', () => {
  let component: TceComponent;
  let fixture: ComponentFixture<TceComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ TceComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(TceComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
